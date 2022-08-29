function [ fr_tmp, fK_tmp, fB_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, x_interp2_tmp,...
    xp_bias_offset_tmp, xr_ideal_tmp, xs_ideal_tmp, xp0_tmp, yp0_tmp, fer_tmp, x0, y0, err_mm, err_deg, mov_dir, ...
    vymax, tvymax, vtmax, tvtmax, Fcom, rt, mtime] = ...
    reject_movements_0924_2021a(fr_tmp, fK_tmp, fB_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, x_interp2_tmp,...
    xp_bias_offset_tmp, xr_ideal_tmp, xs_ideal_tmp, xp0_tmp, yp0_tmp, fer_tmp, x0, y0, err_mm, err_deg, mov_dir, ...
    vymax, tvymax, vtmax, tvtmax, Fcom, rt, mtime,...
    dat, tgt_all, filter_err_flag)

%%%NOTE: When calculating the percentage of trials filtered, I should
%%%account for true differences in the experiment length depending on the
%%%tgt file
%This version also optionally filters trials maxed on error from the ideal
%path (see 0922 or older for versions that didnt do this)
 
% Set threshold parameters
vmax_thresh=0.8;
vmin_thresh=0.2; %try 0.25 and 0.5
mtmax=2;%
pathlenmin = 0.08;
pathlenmax = 0.2;
rtmax = 2 * 1000; %reaction times are passed in msec format
yrms_thresh = 10; %8-10 seems to be a good threshold
minpos_min = -0.0015*1000; %starting point (or minimum) of y-position
minpos_max = 0.0015*1000; %starting point (or maximum) of y-position
maxf_thresh = 18; %robot is not capable of commanding more than 18N total
maxdtheta = 1.5; %maximum deviation from imposed direction on vEC trials
maxvel_delta_thresh = 20; %threshold for the discrepancy between occurrence of trajectory midpoint and max velocity (in samples)
max_err_thresh = 0.8; %maximum allowable error observed on a FC trial
 
% Preallocate
[nmov,nsub] = size(dat.good);
nmov = nmov/2; %forward movements only, so divide by 2
mov_ok = ones(nmov,nsub);

%%%get forward movement data from the dat file
mtime_real = dat.mtime_real(2:2:end,:);
pathlen = dat.pathlen(2:2:end,:);

% Note bad movements
mov_ok(mtime_real>mtmax)=0;
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

%optionally filter trials with large errors from ideal
if filter_err_flag
   fer2 = fer_tmp; 
   fer2(eidx) = {nan(nsamples,nsub)};
   fer = cell2mat(shiftdim(fer2,-2));
   max_err = shiftdim(max(abs(fer)),1)';
   mov_ok(max_err>max_err_thresh)=0;
   max_err_id = nansum(max_err(:)>max_err_thresh(:));
   disp(['Rejected ' num2str(max_err_id/sum(~isnan(max_err(:)))*100) '% of trials with errors above threshold'])
end

% Detect deviation away from intended movement angle on vEC trials, and
% filter based on this deviation
xp = cell2mat(shiftdim(xp_tmp,-2));

%for each trial, find the angle in the middle of the movement
total_displ = sqrt(xp.^2 + yp.^2);

xmid = nan(nmov, nsub);
ymid = nan(nmov, nsub);
tmid = nan(nmov, nsub);
for imov = 1:nmov
    for isub = 1:nsub
        mid_idx = find(squeeze(total_displ(:,isub,imov))>50,1,'first');
        if isempty(mid_idx), mid_idx = ceil(nsamples/2); end
        tmid(imov,isub) = mid_idx;
        xmid(imov,isub) = xp(mid_idx,isub,imov);
        ymid(imov,isub) = yp(mid_idx,isub,imov);
    end
end
% mid_idx = ceil(nsamples/2);
% xmid = shiftdim(xp(mid_idx,:,:),1)' - x0;
% ymid = shiftdim(yp(mid_idx,:,:),1)' - y0;
movang = atand(ymid./xmid);

%for now I wont worry about filtering based on movement angles, but I'll
%need to figure something out later
% movang = 90-movang; %center the angles where striaght ahead is 0 degrees
% movang(movang>90) = movang(movang>90)-180;
% iang = tgt_all(:,end); %imposed angle
% iang(tgt_all(:,4)==0) = nan; %null or FF trial
% iang(tgt_all(:,4)==1) = 0; %imposed error is 0 on standard EC
% iang(tgt_all(:,4)==3) = nan; %interpolated clamp trial
% iang_all = repmat(iang,1,size(movang,2));
% mov_ok(abs(iang_all-movang)>maxdtheta) = 0;
% dtheta_id = sum(abs(iang_all(:)-movang(:))>maxdtheta);
% disp(['dtheta rejected ' num2str(dtheta_id/sum(~isnan(iang_all(:)))*100,2) '% of the vEC movements'])

%filter trials in which the max velocity occurs too long after the midpoint
%of the movement
maxvel_delta = abs(tvtmax-tmid);
mov_ok(maxvel_delta>maxvel_delta_thresh) = 0;
 
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
x0(mov_ok==0) = nan;
y0(mov_ok==0) = nan;
err_mm(mov_ok==0) = nan;
err_deg(mov_ok==0) = nan;
mov_dir(mov_ok==0) = nan;
vymax(mov_ok==0) = nan;
tvymax(mov_ok==0) = nan;
vtmax(mov_ok==0) = nan;
tvtmax(mov_ok==0) = nan;
Fcom(mov_ok==0) = nan;
rt(mov_ok==0) = nan;
mtime(mov_ok==0) = nan;

padnan_sub = nan(nsamples,1);

for k1=1:nmov
    for k2=1:nsub
        if mov_ok(k1,k2)==0
            if ~isempty(fr_tmp{k1})
                fr_tmp{k1}(:,k2) = padnan_sub;
                fK_tmp{k1}(:,k2) = padnan_sub;
                fB_tmp{k1}(:,k2) = padnan_sub;
            end            
            vx_tmp{k1}(:,k2) = padnan_sub;
            vy_tmp{k1}(:,k2) = padnan_sub;
            yp_tmp{k1}(:,k2) = padnan_sub;
            xp_tmp{k1}(:,k2) = padnan_sub;
            ms_tmp{k1}(:,k2) = NaN;
            x_interp2_tmp{k1}(:,k2) = padnan_sub;
            xp_bias_offset_tmp{k1}(:,k2) = padnan_sub;
            xr_ideal_tmp{k1}(:,k2) = padnan_sub;
            xs_ideal_tmp{k1}(:,k2) = padnan_sub;
            xp0_tmp{k1}(:,k2) = padnan_sub;
            yp0_tmp{k1}(:,k2) = padnan_sub;
            fer_tmp{k1}(:,k2) = padnan_sub;
        end
    end
end

end