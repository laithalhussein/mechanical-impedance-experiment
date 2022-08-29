function [ fr_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, err_mm, err_deg, mov_dir, ...
    vymax, vtmax, Fcom, rt, mtime, ideal_err_mm, ideal_err_deg, ideal_max_err] = ...
    reject_movements_0101_2022a(fr_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, err_mm, err_deg, mov_dir, ...
    vymax, vtmax, Fcom, rt, mtime, ideal_err_mm, ideal_err_deg, ideal_max_err,...
    dat, tgt_all)
 
% Set threshold parameters
vmax_thresh=1;
vmin_thresh=0.2; %try 0.25 and 0.5
mtmax=2;%
pathlenmin = 0.1*1000;
pathlenmax = 0.3*1000;
rtmax = 2; %reaction times are passed in msec format
yrms_thresh = 20; %8-10 seems to be a good threshold
minpos_min = -0.0015*1000; %starting point (or minimum) of y-position
minpos_max = 0.0015*1000; %starting point (or maximum) of y-position
maxf_thresh = 18; %robot is not capable of commanding more than 18N total
maxdtheta = 1.5; %maximum deviation from imposed direction on vEC trials
 
% Preallocate
[nmov,nsub] = size(dat.good);
nmov = nmov/2; %forward movements only, so divide by 2
mov_ok = ones(nmov,nsub);

%%%get forward movement data from the dat file
pathlen = dat.pathlen(2:2:end,:);

% Note bad movements
mov_ok(mtime>mtmax)=0;
mov_ok(vtmax>vmax_thresh)=0;
mov_ok(vtmax<vmin_thresh)=0;
mov_ok(rt>rtmax)=0;
mov_ok(pathlen<pathlenmin)=0;
mov_ok(pathlen>pathlenmax)=0;

% Note that in the robot code, the actual restriction is 1.5 mm from the starting
% position.
yp = cell2mat(shiftdim(yp_tmp,-2));
minp = shiftdim(min(yp),1)'; %shiftdim gets rid of singleton, and then transposing puts subjects on columns
mov_ok(minp>minpos_max)=0;
mov_ok(minp<minpos_min)=0;

%similarly detect erroneous forces
%the force cell array is empty when no force data is available, so pad with
%nans to keep indexing consistent
fr2 = fr_tmp;
eidx = cellfun('isempty',fr2); %want to fill trials without force data with NaNs
nsamples = size(yp_tmp{1},1);
fr2(eidx) = {nan(nsamples,nsub)};
fr = cell2mat(shiftdim(fr2,-2));
maxf = shiftdim(max(abs(fr)),1)';
mov_ok(maxf>maxf_thresh)=0;
maxf_id = nansum(maxf(:)>maxf_thresh(:));
disp(['Rejected ' num2str(maxf_id/sum(~isnan(maxf(:)))*100) '% of trials with forces above threshold'])

% Detect deviation away from intended movement angle on vEC trials, and
% filter based on this deviation
xp = cell2mat(shiftdim(xp_tmp,-2));

%for each trial, find the angle in the middle of the movement
total_displ = sqrt(xp.^2 + yp.^2);

xmid = nan(nmov, nsub);
ymid = nan(nmov, nsub);
tmid = nan(nmov, nsub);
%NOTE: I should maybe account for the initial starting position
for imov = 1:nmov
    for isub = 1:nsub
        mid_idx = find(squeeze(total_displ(:,isub,imov))>50,1,'first');
        if isempty(mid_idx), mid_idx = ceil(nsamples/2); end
        tmid(imov,isub) = mid_idx;
        xmid(imov,isub) = xp(mid_idx,isub,imov);
        ymid(imov,isub) = yp(mid_idx,isub,imov);
    end
end
movang_tmp = atand(ymid./xmid);
movang_tmp = 90-movang_tmp; %center the angles where striaght ahead is 0 degrees
movang_tmp(movang_tmp>90) = movang_tmp(movang_tmp>90)-180;
movang2 = squeeze(err_deg(:,:,1)); %all are good?
movang3 = squeeze(mov_dir(:,:,1)); %4.4

movang = movang3;

%for now I wont worry about filtering based on movement angles, but I'll
%need to figure something out later
iang = tgt_all(:,end); %imposed angle
iang(tgt_all(:,7)==0) = nan; %null or FF trial
iang(tgt_all(:,7)==1) = 0; %imposed error is 0 on standard EC
iang_all = repmat(iang,1,size(movang,2));
mov_ok(abs(iang_all-movang)>maxdtheta) = 0;
dtheta_id = sum(abs(iang_all(:)-movang(:))>maxdtheta);
disp(['dtheta rejected ' num2str(dtheta_id/sum(~isnan(iang_all(:)))*100,2) '% of the vEC movements'])

%detect if the movement went backwards any
 
bad_before = sum(mov_ok(:)==0);
% Figure out how far a movement's yp is from that participant's typical
% yp (median of yp after trial 100). Reject movements where the RMS of the error
% exceeds the threshold. Super weird trials get rejected, though its also a
% tiny fraction
yp_standard = nanmedian(yp(:,:,101:end),3);
yp_standard_rep = repmat(yp_standard,1,1,nmov);
% yse = bsxfun(@(x,y) ((x-y).^2)/sum(~isnan(x)&~isnan(y)), yp,yp_standard); %normalized squared error
% yrms = shiftdim(sqrt(nansum(yse)),1)';
yrms = shiftdim(rms(yp-yp_standard_rep),1)';
mov_ok(yrms>yrms_thresh)=0;
 
bad_after = sum(mov_ok(:)==0);
rms_detected = bad_after-bad_before;
ntot = length(mov_ok(:));
identified = sum(yrms(:)>yrms_thresh);
disp(['Distance from standard position profile flagged ' num2str(identified/ntot*100,2) '% of movements'])

%now we are ready to actually remove those trials (nan them!)
total_accepted = 1 - sum(mov_ok(:)==0)/length(mov_ok(:));
disp(['We accepted ', num2str(total_accepted*100), '% of trials']);

%first deal with 3D matrices
err_mm1_tmp = squeeze(err_mm(:,:,1));
err_mm2_tmp = squeeze(err_mm(:,:,2));
err_deg1_tmp = squeeze(err_deg(:,:,1));
err_deg2_tmp = squeeze(err_deg(:,:,2));
mov_dir1_tmp = squeeze(mov_dir(:,:,1));
mov_dir2_tmp = squeeze(mov_dir(:,:,2));
Fcom1_tmp = squeeze(Fcom(:,:,1));
Fcom2_tmp = squeeze(Fcom(:,:,2));

err_mm1_tmp(mov_ok==0) = nan;
err_mm2_tmp(mov_ok==0) = nan;
err_deg1_tmp(mov_ok==0) = nan;
err_deg2_tmp(mov_ok==0) = nan;
mov_dir1_tmp(mov_ok==0) = nan;
mov_dir2_tmp(mov_ok==0) = nan;
Fcom1_tmp(mov_ok==0) = nan;
Fcom2_tmp(mov_ok==0) = nan;

err_mm(:,:,1) = err_mm1_tmp;
err_mm(:,:,2) = err_mm2_tmp;
err_deg(:,:,1) = err_deg1_tmp;
err_deg(:,:,2) = err_deg2_tmp;
mov_dir(:,:,1) = mov_dir1_tmp;
mov_dir(:,:,2) = mov_dir2_tmp;
Fcom(:,:,1) = Fcom1_tmp;
Fcom(:,:,2) = Fcom2_tmp;

%then 2D matrices
vymax(mov_ok==0) = nan;
vtmax(mov_ok==0) = nan;
rt(mov_ok==0) = nan;
mtime(mov_ok==0) = nan;

%then cell arrays
padnan_sub = nan(nsamples,1);
for k1=1:nmov
    for k2=1:nsub
        if mov_ok(k1,k2)==0
            if ~isempty(fr_tmp{k1})
                fr_tmp{k1}(:,k2) = padnan_sub;
            end            
            vx_tmp{k1}(:,k2) = padnan_sub;
            vy_tmp{k1}(:,k2) = padnan_sub;
            yp_tmp{k1}(:,k2) = padnan_sub;
            xp_tmp{k1}(:,k2) = padnan_sub;
            ms_tmp{k1}(:,k2) = padnan_sub;
        end
    end
end

return