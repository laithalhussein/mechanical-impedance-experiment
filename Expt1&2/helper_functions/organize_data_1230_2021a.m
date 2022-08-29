function [dat,info]=organize_data_1230_2021a(draw, iraw)
%This function should (1) align the data, (2) extract the indices we want, (3) filter the data, (4) return the ILF for each subject

close all;
num_subjects = length(iraw.sublist);
%num_samples = iraw.num_samples;

filter_trials_flag = 1;

info.FFMAG = iraw.FFMAG;
info.sublist = iraw.sublist;
EC_err_size = 5; %in degrees

max_trials = size(draw.fyr,2);

%pre-allocate array to hold the index of the block breaks for each participant
num_blocks = length(iraw.blocksize);
new_block_dummy = nan(1,num_blocks);

%% get the indices of interest
% fam_block = 1;
% baseline_block = [2,3];
% training1_block = [4,5];
% test1_block = [6,7,8,9,10,11];
% test2_block = [12,13,14,15,16,17];
% training3_block = [18,19];

fam_idx = [1:60];
baseline_idx = [fam_idx(end)+ 1 : fam_idx(end) + 200];
training1_idx = [baseline_idx(end) + 1 : baseline_idx(end) + 300];
test1_idx = [training1_idx(end) + 1 : training1_idx(end) +  990];
test2_idx = [test1_idx(end) + 1 : test1_idx(end) + 990];
training2_idx = [test2_idx(end) + 1: test2_idx(end) + 300];

%concentrate on forward trials...
fam_idx = fam_idx(2:2:end)/2;
baseline_idx = baseline_idx(2:2:end)/2;
training1_idx = training1_idx(2:2:end)/2;
test1_idx = test1_idx(2:2:end)/2;
test2_idx = test2_idx(2:2:end)/2;
training2_idx = training2_idx(2:2:end)/2;

%save experiment sequence information
info.exp_seq.fam_idx = fam_idx;
info.exp_seq.baseline_idx = baseline_idx;
info.exp_seq.training1_idx = training1_idx;
info.exp_seq.test1_idx = test1_idx;
info.exp_seq.test2_idx = test2_idx;
info.exp_seq.training2_idx = training2_idx;

info.tpb = [length(fam_idx), length(baseline_idx), length(training1_idx), length(test1_idx), length(test2_idx),...
    length(training2_idx)];

tgt_all = cell(num_subjects,1);
for k1=1:num_subjects
    
    idx.baseline.all(:,k1) = baseline_idx;
    idx.training1(:,k1) = training1_idx;
    idx.test1.all(:,k1) = test1_idx;
    idx.training2(:,k1) = training2_idx;
    idx.test2.all(:,k1) = test2_idx;
    
    sub_tgt_tmp = iraw.tgt_dat(k1,:);
    sub_tgt_all = vertcat(sub_tgt_tmp{:});
    sub_tgt_all_fwrd = sub_tgt_all(2:2:end,:);
    tgt_all{k1} = sub_tgt_all_fwrd;
    
    %get all EC trials
    ec_tmp1 = find(sub_tgt_all_fwrd(:,7) == 1);
    
    %get all vEC trials
    ec_tmp2 = find(sub_tgt_all_fwrd(:,7) == 2 ); %not counting the vEC triplet
    
    %get EC trials during baseline
    idx.baseline.ec(:,k1) = ec_tmp1(ec_tmp1 <= idx.baseline.all(end,k1));
    
    %isolate vEC trials during the test periods
    idx.vEC.test1(:,k1) = ec_tmp2( ec_tmp2<=idx.test1.all(end,k1) & ec_tmp2>=idx.test1.all(1,k1));
    idx.vEC.test2(:,k1) = ec_tmp2( ec_tmp2<=idx.test2.all(end,k1) & ec_tmp2>=idx.test2.all(1,k1));
    
    %find the triplets for each case
    FF_triplets_all_tmp = find(sub_tgt_all_fwrd(:,2)~=0); %FF
    FF_triplets_P_tmp = find(sub_tgt_all_fwrd(:,2)==-1);
    FF_triplets_N_tmp = find(sub_tgt_all_fwrd(:,2)==1);
    
    FC_triplets_all_tmp = find(abs(sub_tgt_all_fwrd(:,end))==EC_err_size); %FC
    FC_triplets_P_tmp = find(sub_tgt_all_fwrd(:,end)==EC_err_size);
    FC_triplets_N_tmp = find(sub_tgt_all_fwrd(:,end)==-EC_err_size);
    
    idx.pert.FF.all(:,k1) = FF_triplets_all_tmp;
    idx.pert.FC.all(:,k1) = FC_triplets_all_tmp;
    
    idx.pert.FF.P(:,k1) = FF_triplets_P_tmp;
    idx.pert.FC.P(:,k1) = FC_triplets_P_tmp;
    
    idx.pert.FF.N(:,k1) = FF_triplets_N_tmp;
    idx.pert.FC.N(:,k1) = FC_triplets_N_tmp;
    
    idx.pre.FF.all(:,k1) = FF_triplets_all_tmp - 1;     idx.post.FF.all(:,k1) = FF_triplets_all_tmp + 1;
    idx.pre.FC.all(:,k1) = FC_triplets_all_tmp - 1;     idx.post.FC.all(:,k1) = FC_triplets_all_tmp + 1;
    
    idx.pre.FF.P(:,k1) = FF_triplets_P_tmp - 1;     idx.post.FF.P(:,k1) = FF_triplets_P_tmp + 1;
    idx.pre.FC.P(:,k1) = FC_triplets_P_tmp - 1;     idx.post.FC.P(:,k1) = FC_triplets_P_tmp + 1;
    
    idx.pre.FF.N(:,k1) = FF_triplets_N_tmp - 1;     idx.post.FF.N(:,k1) = FF_triplets_N_tmp + 1;
    idx.pre.FC.N(:,k1) = FC_triplets_N_tmp - 1;     idx.post.FC.N(:,k1) = FC_triplets_N_tmp + 1;
    
    %find the sign of the applied errors during the triplets
    %IMPORTANT: Positive ECs induce leftward displacements, whereas a +FF induces rightward displacements    
    err_sign.FF(:,k1) = sub_tgt_all_fwrd(FF_triplets_all_tmp,2);
    err_sign.EC(:,k1) = sign(sub_tgt_all_fwrd(FC_triplets_all_tmp,end));
    
    %combine all vEC trials together
    idx.vEC.all(:,k1) = ec_tmp2;
    
    %get all non-baseline EC trials
    idx.EC.all(:,k1) = ec_tmp1(ec_tmp1 > idx.baseline.all(end,k1));    
end

info.idx = idx;
info.tgt_all = tgt_all;
info.err_sign = err_sign;

%% align the data
%pre-allocate cell arrays to hold important forward trial data
[vx_tmp, vy_tmp, yp_tmp, xp_tmp, fr_tmp, ms_tmp] = deal(cell(max_trials/2, 1));
%[vxa_tmp, vya_tmp, ypa_tmp, xpa_tmp, fra_tmp, msa_tmp] = deal(cell(max_trials/2, 1));

%save the kinematic error or movement direction at: midpoint of movement and endpoint of matching region
err_mm = nan(max_trials/2, num_subjects, 2); %last dimension indexes the point into the movement
err_deg = nan(max_trials/2, num_subjects, 2); %save the error in degrees as well
mov_dir = nan(max_trials/2, num_subjects, 2);

ideal_err_mm = nan(max_trials/2, num_subjects, 2); %error based on the ideal path
ideal_err_deg = nan(max_trials/2, num_subjects, 2);
ideal_max_err = nan(max_trials/2, num_subjects);

%save the commanded force at given points into the movement
Fcom = nan(max_trials/2, num_subjects, 2);

%save the max velocity data
vymax = draw.vxmax(2:2:end,:);
vtmax = draw.vtmax(2:2:end,:);

%save movement time data
mtime = draw.mtime_real(2:2:end,:); %put it in ms
rt = draw.rt(2:2:end,:);

tgt_win = [-79:0];
mo_win = [-45:35];
end_win = [0:200];
win = tgt_win;

short_movements = 0;
for k2=1:num_subjects
    cc=1;
    for k3=2:2:max_trials
        tgt_idx = find(draw.pxr{k3}(:,k2) > 0.1, 1, 'first') - 1;
        pillow_idx = find(draw.pxr{k3}(:,k2) > 0.11, 1, 'first') - 1;
        mov_onset_idx = find(draw.ms{k3}(:,k2)==4,1,'first');
        
        total_displacement = sqrt(draw.pxr{k3}(:,k2).^2 + draw.pyr{k3}(:,k2).^2);
        midpoint_idx = find(total_displacement>0.05, 1, 'first') - 1;
        endpoint_idx = find(total_displacement>0.1, 1, 'first') - 1;
        
        force_tmp = draw.fyr{k3};
        
        %we can have cases in the early trials where there is no force data saved (these are just null trials anyway)
        if size(force_tmp,2)~=num_subjects, force_tmp(:,end:num_subjects) = NaN; end
        
        %sometimes the movement can be short, and we dont have enough data
        %before the pillow to subsample the movement
        %         if (~isempty(tgt_idx) & tgt_idx<length(win)) | (tgt_idx + length(end_win) > size(force_tmp,1)) | mov_onset_idx<=abs(mo_win(1))
        if draw.good(k3,k2)==0 | isempty(pillow_idx) | mov_onset_idx<=abs(mo_win(1)) |...
                (tgt_idx + length(end_win) > size(force_tmp,1)) | endpoint_idx<=midpoint_idx
            vx_tmp{cc}(:,k2) = nan+win;
            vy_tmp{cc}(:,k2) = nan+win;
            yp_tmp{cc}(:,k2) = nan+win;
            xp_tmp{cc}(:,k2) = nan+win;
            fr_tmp{cc}(:,k2) = nan+win;
            ms_tmp{cc}(:,k2) = nan+win;
            if ~isempty(force_tmp)
                fr_tmp{cc}(:,k2) = nan+win;
            else
                fr_tmp{cc} = [];
            end
            short_movements = short_movements+1;            
        else
            %if isempty(tgt_idx), tgt_idx = abs(min(win)) + 1; end
            
            vy_tmp{cc}(:,k2) = draw.vxr{k3}(tgt_idx+win,k2);
            vx_tmp{cc}(:,k2) = draw.vyr{k3}(tgt_idx+win,k2);
            yp_tmp{cc}(:,k2) = draw.pxr{k3}(tgt_idx+win,k2) * 1000;
            xp_tmp{cc}(:,k2) = draw.pyr{k3}(tgt_idx+win,k2) * 1000;
            ms_tmp{cc}(:,k2) = draw.ms{k3}(tgt_idx+win,k2);
            
            %find the initial start position
            start_idx = find(draw.ms{k3}(:,k2)==2,1,'first');
            if isempty(start_idx), 
                start_idx = find(draw.ms{k3}(:,k2)==3,1,'first');
            end
            x0 = draw.pyr{k3}(start_idx,k2) * 1000;
            y0 = draw.pxr{k3}(start_idx,k2) * 1000;
            
            %calculate error based on ideal movement
            y_len = 100-y0; %calculate the length of the movement to the end target (longitudinally)
            y_pos = yp_tmp{cc}(:,k2)-y0; %get the actual y position after centering
            
            %normalize y-pos by length of the movement,
            %and then subtract from 1 to get a function that describes how y *should* change over time
            y_frac = 1-y_pos/y_len;
            
            %now that we know how y should ideally change, we can see how x
            %should ideally change by multiplying the above function to the starting x position
            x_ideal = y_frac*x0;% the ideal x path
            
            %calculate the error
            x_err = xp_tmp{cc}(:,k2)-x_ideal; %y-deviation from the straight line
            
            %recalculate some indices based on aligned data
            td_tmp = sqrt(xp_tmp{cc}(:,k2).^2 + yp_tmp{cc}(:,k2).^2);
            tmid_tmp = find(td_tmp>50,1,'first')-1;
            mo_tmp = find(ms_tmp{cc}(:,k2)==4,1,'first')-1;
            
            %get raw errors/mov direction
            err_mm(cc,k2,1) = draw.pyr{k3}(midpoint_idx,k2) * 1000;            
            err_mm(cc,k2,2) = draw.pyr{k3}(tgt_idx,k2) * 1000;
                        
            if ~isempty(mo_tmp) & mo_tmp~=0 & mo_tmp<tmid_tmp & tmid_tmp~=0
                %center the angles where striaght ahead is 0 degrees
                e1_tmp = 90 - atan2d(y_pos(tmid_tmp)-y_pos(mo_tmp), xp_tmp{cc}(tmid_tmp,k2)-xp_tmp{cc}(mo_tmp,k2));
                e2_tmp = 90 - atan2d(y_pos(end)-y_pos(mo_tmp), xp_tmp{cc}(end,k2)-xp_tmp{cc}(mo_tmp,k2));
                if e1_tmp>90, e1_tmp = e1_tmp-180; end
                if e2_tmp>90, e2_tmp = e2_tmp-180; end
                mov_dir(cc,k2,1) = e1_tmp;
                mov_dir(cc,k2,2) = e2_tmp;
                
                err_deg(cc,k2,1) = atand(x_err(tmid_tmp)/y_pos(tmid_tmp));
                err_deg(cc,k2,2) = atand(x_err(end)/y_pos(end));
                
                %save ideal errors
                ideal_err_mm(cc,k2,1) = x_err(tmid_tmp);
                ideal_err_mm(cc,k2,2) = x_err(end);
                ideal_err_deg(cc,k2,1) = atand(x_err(tmid_tmp)/y_pos(tmid_tmp));
                ideal_err_deg(cc,k2,2) = atand(x_err(end)/y_pos(end));
                ideal_max_err(cc,k2) = max(x_err);
            end
            if isempty(force_tmp)
                fr_tmp{cc}= [];
            else
                fr_tmp{cc}(:,k2) = force_tmp(win+tgt_idx, k2);
                if ~isempty(tmid_tmp) && tmid_tmp~=0, Fcom(cc,k2,1) = fr_tmp{cc}(tmid_tmp,k2); end
                Fcom(cc,k2,2) = fr_tmp{cc}(end,k2);
            end
        end
        cc=cc+1;
    end
end

%% next filter the data
%%%pass the data through our filtering function, but AFTER data is
%%%aligned...?
if filter_trials_flag
    [fr_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, err_mm,...
        err_deg, mov_dir, vymax, vtmax, Fcom, rt, mtime,...
        ideal_err_mm, ideal_err_deg, ideal_max_err] = ...
        reject_movements_0101_2022a(fr_tmp, vx_tmp, vy_tmp, yp_tmp, xp_tmp, ms_tmp, err_mm,...
        err_deg, mov_dir, vymax, vtmax, Fcom, rt, mtime, ideal_err_mm, ideal_err_deg, ideal_max_err,...
        draw, sub_tgt_all_fwrd);
end

%% calculate ILF
ILF_norm = nan(max_trials/2, num_subjects);
ILF_raw = nan(max_trials/2, num_subjects);
F_ILF = @mean;
for k1=1:max_trials/2
    if ~isempty(fr_tmp{k1})
        for k2=1:num_subjects
            ILF_norm(k1,k2) = F_ILF(fr_tmp{k1}(:,k2)) / F_ILF( vy_tmp{k1}(:,k2) * info.FFMAG);
            %if ~all(isnan(fr_tmp{k1}(:,k2))) & ~all(fr_tmp{k1}(:,k2)==fr_tmp{k1}(1,k2)) & k1>200, keyboard; end
            ILF_raw(k1,k2) = F_ILF(fr_tmp{k1}(:,k2));
        end
    end
end

%% save data for dat structure, and convert cell data to 3D matrix

%switch 3rd and 1st dimension to put data in trials x subject x sample format
dat.vx_tmp = permute(cell2mat(shiftdim(vx_tmp,-2)),[3,2,1]);
dat.vy_tmp = permute(cell2mat(shiftdim(vy_tmp,-2)),[3,2,1]);
dat.yp_tmp = permute(cell2mat(shiftdim(yp_tmp,-2)),[3,2,1]);
dat.xp_tmp = permute(cell2mat(shiftdim(xp_tmp,-2)),[3,2,1]);
dat.fr_tmp = permute(cell2mat(shiftdim(fr_tmp,-2)),[3,2,1]);
dat.ms_tmp = permute(cell2mat(shiftdim(ms_tmp,-2)),[3,2,1]);

dat.vymax = vymax;
dat.vtmax = vtmax;
dat.err_mm = err_mm;
dat.err_deg = err_deg;
dat.mov_dir = mov_dir;
dat.ideal_err_mm = ideal_err_mm;
dat.ideal_err_deg = ideal_err_deg;
dat.ideal_max_err = ideal_max_err;
dat.Fcom = Fcom;
dat.rt = rt;
dat.mtime = mtime;
dat.ILF_norm = ILF_norm;
dat.ILF_raw = ILF_raw;

info.num_samples = length(win);
info.num_original_samples = iraw.num_samples;

return


