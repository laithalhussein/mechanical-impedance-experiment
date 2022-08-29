function [dat,info]=organize_data_1205_2021a(draw, iraw)
%This function should (1) align the data, (2) extract the indices we want, (3) filter the data, (4) return the ILF for each subject

close all;
num_subjects = length(iraw.sublist);
num_samples = iraw.num_samples;

filter_trials_flag = 1;

%define a flag to optionally filter the force data with low-pass filtering
filt_force_flag = 0;

num_triplets = 33; %there are 48 triplets for each test phase
EC_err_size = 5; %in degrees
EC_err_mm = tand(EC_err_size)*100;

%begin to define dat and info
dat = [];
info.FFMAG = 4;
info.sublist = iraw.sublist;

max_trials = size(draw.fyr,2);

%pre-allocate cell arrays to hold important forward trial data
[vx_tmp, vy_tmp, yp_tmp, xp_tmp, fr_tmp, ms_tmp] = deal(cell(max_trials/2, 1));
[vxa_tmp, vya_tmp, ypa_tmp, xpa_tmp, fra_tmp, msa_tmp] = deal(cell(max_trials/2, 1));

%save FF size
FF_size = nan(1,num_subjects);

%save the kinematic error or movement direction at: midpoint of movement, endpoint of matching region, time of max velocity
err_mm = nan(max_trials/2, num_subjects, 3); %last dimension indexes the point into the movement
err_deg = nan(max_trials/2, num_subjects, 3); %save the error in degrees as well
mov_dir = nan(max_trials/2, num_subjects, 3);

%save the imposed error
ierr_mm = nan(max_trials/2);

%save the commanded force at given points into the movement (will be used to estimate stiffness later)
Fcom = nan(max_trials/2, num_subjects, 3);
ILF = nan(max_trials/2, num_subjects);

%save the max velocity data
vymax = draw.vxmax(2:2:end,:);
vtmax = draw.vtmax(2:2:end,:);

%save movement time data
%mtime = draw.mtime_real(2:2:end,:);
mtime = draw.mtime(2:2:end,:)*0.005;

%pre-allocate arrays to save the time at which max velocity occurs
tvymax = nan(max_trials/2, num_subjects);
tvtmax = nan(max_trials/2, num_subjects);

%define flags that determine if loop through a given trial
num_bad_pathlen = 0;
num_bad_start = 0;

%save reaction times
%rt = draw.rt(2:2:end,:)*1000; %in ms
rt = mtime*NaN;

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
    
    idx.baseline.all(k1,:) = baseline_idx;    
    idx.training1(k1,:) = training1_idx;
    idx.test1.all(k1,:) = test1_idx;    
    idx.training2(k1,:) = training2_idx;
    idx.test2.all(k1,:) = test2_idx;
    
    sub_tgt_tmp = info.tgt_dat(k1,:);
    sub_tgt_all = vertcat(sub_tgt_tmp{:});
    sub_tgt_all_fwrd = sub_tgt_all(2:2:end,:);
    tgt_all{k1} = sub_tgt_all_fwrd;
    
    %get all EC trials
    ec_tmp1 = find(sub_tgt_all_fwrd(:,7) == 1);
    
    %get all vEC trials
    ec_tmp2 = find(sub_tgt_all_fwrd(:,7) == 2 ); %not counting the vEC triplet
    
    %get EC trials during baseline
    idx.baseline.ec(k1,:) = ec_tmp1(ec_tmp1 <= idx.baseline.all(k1,end));
    
    %isolate vEC trials during the test periods
    idx.vEC.test1(k1,:) = ec_tmp2( ec_tmp2<=idx.test1.all(k1,end) & ec_tmp2>=idx.test1.all(k1,1) );
    idx.vEC.test2(k1,:) = ec_tmp2( ec_tmp2<=idx.test2.all(k1,end) & ec_tmp2>=idx.test2.all(k1,1) );
    
    %find the triplets for each case
    triplets_tmp1 = find(sub_tgt_all_fwrd(:,2)~=0); %FF
    triplets_tmp2 = find(abs(sub_tgt_all_fwrd(:,end))==EC_err_size); %EC
    
    idx.err_FF.all = triplets_tmp1;
    idx.err_EC.all = triplets_tmp2;
    
    idx.pre.err_FF(k1,:) = triplets_tmp1 -1;     idx.post.err_FF(k1,:) = triplets_tmp1+1;
    idx.pre.err_EC(k1,:) = triplets_tmp2 -1;     idx.post.err_EC(k1,:) = triplets_tmp2 +1;
    
    %find the sign of the applied errors during the triplets
    %IMPORTANT: Positive ECs induce leftward displacements, whereas a +FF induces rightward displacements    
    err_sign.FF(k1,:) = sub_tgt_all_fwrd(triplets_tmp1,2);
    err_sign.EC(k1,:) = sign(sub_tgt_all_fwrd(triplets_tmp2,end));
    
    %combine all vEC trials together
    idx.vEC.all(k1,:) = ec_tmp2;
    
    %get all non-baseline EC trials
    idx.EC.all(k1,:) = ec_tmp1(ec_tmp1 > idx.baseline.all(k1,end));    
end

%% align the data
cc=1;
for ntrial=2:2:max_trials
    
    pathlen_flag = 0;
    if end_mov_align==0
        if ~all(isnan(draw.ms{ntrial}(:,nsub))) && max(draw.pxr{ntrial}(:,nsub)*1000) < end_mov_thresh,
            num_bad_pathlen = num_bad_pathlen+1;
            pathlen_flag = 1;
        end
    end
    
    bad_start = 0;
    if end_mov_align==0
        if all(draw.pxr{ntrial}(1,nsub)*1000 > end_mov_thresh) %if the first x samples are outside bounds we would reasonably expect, exclude the trial
            bad_start = 1;
            num_bad_start = num_bad_start+1;
            %keyboard;
        end
    end
    
    if ~all(isnan(draw.ms{ntrial}(:,nsub))) && ~pathlen_flag && ~bad_start %dont bother if its a bad trial...
        
        %find movement onset and index of tgt appearance
        %         mov_onset_idx = find(dat.ms{ntrial}(:,nsub)==4,1,'first');
        %         tgt_show_idx = find(dat.ms{ntrial}(:,nsub)==2,1,'last');
        
        %save force data now for filtering later
        force_tmp = draw.fyr{ntrial};
        
        %save data (and make pos units mm)
        vx_tmp{cc}(:,nsub) = draw.vyr{ntrial}(:,nsub);
        vy_tmp{cc}(:,nsub) = draw.vxr{ntrial}(:,nsub);
        yp_tmp{cc}(:,nsub) = draw.pxr{ntrial}(:,nsub)*1000;
        xp_tmp{cc}(:,nsub) = draw.pyr{ntrial}(:,nsub)*1000;
        ms_tmp{cc}(:,nsub) = draw.ms{ntrial}(:,nsub);
        x_interp2_tmp{cc}(:,nsub) = draw.x_interp2{ntrial}(:,nsub)*1000;
        xp_bias_offset_tmp{cc}(:,nsub) = draw.xp_bias_offset{ntrial}(:,nsub)*1000;
        xp0_tmp{cc}(:,nsub) = draw.xp0{ntrial}(:,nsub)*1000; %note xp0 and yp0 are named correctly
        yp0_tmp{cc}(:,nsub) = draw.yp0{ntrial}(:,nsub)*1000;
        xr_ideal_tmp{cc}(:,nsub) = draw.xr_ideal{ntrial}(:,nsub)*1000;
        xs_ideal_tmp{cc}(:,nsub) = draw.xs_ideal{ntrial}(:,nsub)*1000;
        fer_tmp{cc}(:,nsub) = draw.fer{ntrial}(:,nsub)*1000;
        
        %determine time at which max velocity occurs
        [~,tvymax(cc,nsub)]=max(vy_tmp{cc}(:,nsub));
        vt = sqrt( (vy_tmp{cc}(:,nsub)).^2 + (vx_tmp{cc}(:,nsub)).^2 );
        [~,tvtmax(cc,nsub)]=max(vt);
        
        %%%lets calculate the deviation from the ideal straight line
        %%%between the start and end target throughout different points in
        %%%the movement
        y_len = 100-y0(cc,nsub); %calculate the length of the movement to the end target (longitudinally)
        y_pos = yp_tmp{cc}(:,nsub)-y0(cc,nsub); %get the actual y position after centering
        
        %normalize y-pos by length of the movement,
        %and then subtract from 1 to get a function that describes how y *should* change over time
        y_frac = 1-y_pos/y_len;
        
        %now that we know how y should ideally change, we can see how x
        %should ideally change by multiplying the above function to the starting x position
        x_ideal = y_frac*x0(cc,nsub);% the ideal x path
        
        %calculate the error
        x_err = xp_tmp{cc}(:,nsub)-x_ideal; %y-deviation from the straight line
        
        %determine when the position exceeded the thresholds we're
        %interested in
        total_displacement = sqrt(yp_tmp{cc}(:,nsub).^2 + xp_tmp{cc}(:,nsub).^2);
        mid_idx = find(total_displacement>50,1,'first');
        if end_mov_align==0
            %safer to take the sample before the interpolation ends in case the last sample was corrupted
            yend_idx = find(yp_tmp{cc}(:,nsub)>end_mov_thresh,1,'first')-2;
        else
            yend_idx = num_samples-1;
        end
        
        %save the errors at points of interest
        try
            err_mm(cc,nsub,1) = x_err(mid_idx); %keep in mm
            err_mm(cc,nsub,2) = x_err(yend_idx);
        catch
            keyboard;
        end
        
        %if cc>118, keyboard; end
        %get errors in degrees as well (but use centered y position)
        err_deg(cc,nsub,1) = atand(x_err(mid_idx)/y_pos(mid_idx));
        err_deg(cc,nsub,2) = atand(x_err(yend_idx)/y_pos(yend_idx));
        err_deg(cc,nsub,3) = atand(x_err(tvtmax(cc,nsub))/y_pos(tvtmax(cc,nsub)));
        
        %also calculate movement direction at midpoint and endpoint
        %but with respect to movement onset (then compare to above!)
        mov_onset_idx = find(ms_tmp{cc}(:,nsub)==4,1,'first');
        if ~isempty(mov_onset_idx) & mov_onset_idx<mid_idx
            %center the angles where striaght ahead is 0 degrees
            e1_tmp = 90 - atan2d(y_pos(mid_idx)-y_pos(mov_onset_idx), xp_tmp{cc}(mid_idx,nsub)-xp_tmp{cc}(mov_onset_idx,nsub));
            e2_tmp = 90 - atan2d(y_pos(yend_idx)-y_pos(mov_onset_idx), xp_tmp{cc}(yend_idx,nsub)-xp_tmp{cc}(mov_onset_idx,nsub));
            e3_tmp = 90 - atan2d(y_pos(tvtmax(cc,nsub))-y_pos(mov_onset_idx), xp_tmp{cc}(tvtmax(cc,nsub),nsub)-xp_tmp{cc}(mov_onset_idx,nsub));
            if e1_tmp>90, e1_tmp = e1_tmp-180; end
            if e2_tmp>90, e2_tmp = e2_tmp-180; end
            if e3_tmp>90, e3_tmp = e3_tmp-180; end
            mov_dir(cc,nsub,1) = e1_tmp;
            mov_dir(cc,nsub,2) = e2_tmp;
            mov_dir(cc,nsub,3) = e3_tmp;
        end
        
        if isempty(force_tmp)
            fr_tmp{cc}= [];
            fK_tmp{cc} = [];
            fB_tmp{cc} = [];
        else
            try
                fr_tmp{cc}(:,nsub) = force_tmp(:,nsub);
            catch
                keyboard
            end
            fK_tmp{cc}(:,nsub) = draw.fKr{ntrial}(:,nsub);
            fB_tmp{cc}(:,nsub) = draw.fBr{ntrial}(:,nsub);
            
            Fcom(cc,nsub,1) = fr_tmp{cc}(mid_idx,nsub);
            Fcom(cc,nsub,2) = fr_tmp{cc}(yend_idx,nsub);
            Fcom(cc,nsub,3) = fr_tmp{cc}(tvtmax(cc,nsub),nsub);
        end
    else %if its a bad trial, keep data as nans...
        force_tmp = draw.fyr{ntrial};
        
        if ~isempty(force_tmp)
            fr_tmp{cc}(:,nsub) = force_tmp(:,nsub);
            fK_tmp{cc}(:,nsub) = draw.fKr{ntrial}(:,nsub);
            fB_tmp{cc}(:,nsub) = draw.fBr{ntrial}(:,nsub);
        end
        vx_tmp{cc}(:,nsub) = draw.vyr{ntrial}(:,nsub);
        vy_tmp{cc}(:,nsub) = draw.vxr{ntrial}(:,nsub);
        yp_tmp{cc}(:,nsub) = draw.pxr{ntrial}(:,nsub);
        xp_tmp{cc}(:,nsub) = draw.pyr{ntrial}(:,nsub);
        ms_tmp{cc}(:,nsub) = draw.ms{ntrial}(:,nsub);
        x_interp2_tmp{cc}(:,nsub) = draw.x_interp2{ntrial}(:,nsub);
        xp_bias_offset_tmp{cc}(:,nsub) = draw.xp_bias_offset{ntrial}(:,nsub);
        xp0_tmp{cc}(:,nsub) = draw.xp0{ntrial}(:,nsub);
        yp0_tmp{cc}(:,nsub) = draw.yp0{ntrial}(:,nsub);
        xr_ideal_tmp{cc}(:,nsub) = draw.xr_ideal{ntrial}(:,nsub);
        xs_ideal_tmp{cc}(:,nsub) = draw.xs_ideal{ntrial}(:,nsub);
        fer_tmp{cc}(:,nsub) = draw.fer{ntrial}(:,nsub);
    end
    cc=cc+1;
end

%% next filter the data
%%%pass the data through our filtering function, but AFTER data is
%%%aligned...?
if filter_trials_flag
    [fra_tmp, fKa_tmp, fBa_tmp, vxa_tmp, vya_tmp, ypa_tmp, xpa_tmp, msa_tmp, x_interp2a_tmp,...
        xp_bias_offseta_tmp, xr_ideala_tmp, xs_ideala_tmp, xp0a_tmp, yp0a_tmp, fera_tmp, x0, y0, err_mm,...
        err_deg, mov_dir, vymax, tvymax, vtmax, tvtmax, Fcom, rt, mtime] = ...
        reject_movements_0924_2021a(fra_tmp, fKa_tmp, fBa_tmp, vxa_tmp, vya_tmp, ypa_tmp,...
        xpa_tmp, msa_tmp, x_interp2a_tmp, xp_bias_offseta_tmp, xr_ideala_tmp, xs_ideala_tmp,...
        xp0a_tmp, yp0a_tmp, fera_tmp, x0, y0, err_mm, err_deg, mov_dir, vymax, tvymax, vtmax, tvtmax, Fcom, rt,...
        mtime, draw, ics.tgt_all, filter_err_flag);
end

%% finally we must determine all the indices of interest
for nsub=1:num_subjects
    tgt_tmp = ics.tgt_all{nsub};
    
    idx.fam(:,nsub) = ics.exp_seq{nsub}.fam;
    idx.baseline1(:,nsub) = ics.exp_seq{nsub}.baseline1;
    idx.baseline2(:,nsub) = ics.exp_seq{nsub}.baseline2;
    idx.baseline3(:,nsub) = ics.exp_seq{nsub}.baseline3;
    idx.FF_training(:,nsub) = ics.exp_seq{nsub}.FF_training;
    idx.FF_test(:,nsub) = ics.exp_seq{nsub}.FF_test;
    idx.FC_training(:,nsub) = ics.exp_seq{nsub}.FC_training;
    idx.FC_test(:,nsub) = ics.exp_seq{nsub}.FC_test; %FF and int FC
    
    %find all standard EC trials
    ec_tmp1 = find(tgt_tmp(:,4) == 1);
    
    %find all vEC trials
    ec_tmp2 = find(tgt_tmp(:,4) == 2);
    
    %find baseline EC trials
    idx.baseline_ec(:,nsub) = ec_tmp1(ec_tmp1 <= idx.baseline2(end,nsub));
    
    %get baseline null trials after familiarization
    bln_null_tmp1 = find(tgt_tmp(:,4) == 0 & tgt_tmp(:,2) == 0);
    idx.baseline_null(:,nsub) = bln_null_tmp1(bln_null_tmp1 > idx.fam(end,nsub) & bln_null_tmp1 <= idx.baseline2(end,nsub));
    
    %find the interpolated FC tripelts
    FC_probe_tmp = find(tgt_tmp(:,4) == 3);
    idx.pert.FC.all(:,nsub) = FC_probe_tmp; %interpolated FC
    
    %Separate triplets based on direction
    %Note that positive clamps induce leftward displacements, but positive FFs induce rightward displacements
    
    %%%% first handle FF case
    %NOTE: Might want to update the blow to ensure FF training trials
    %arent captured by mistake
    triplets_tmp_FF = find(abs(tgt_tmp(:,2))==1);
    triplets_tmp_FF_P = find(tgt_tmp(:,2)==1);
    triplets_tmp_FF_N = find(tgt_tmp(:,2)==-1);
    
    %save FF indices for phase 1 of the experiment
    idx.pert.FF.all(:,nsub) = triplets_tmp_FF;
    idx.pert.FF.P(:,nsub) = triplets_tmp_FF_P;
    idx.pert.FF.N(:,nsub) = triplets_tmp_FF_N;
    
    %determine the indices of the refresher trials that surround the FF & FC triplets
    refr_tmp = find((tgt_tmp(:,2)==0 & tgt_tmp(:,4)==0) | (tgt_tmp(:,4)==2)); %refresh trials can only be null or vEC trials
    idx.refr.FF.all(:,nsub) = refr_tmp(refr_tmp>=idx.FF_test(1,nsub) & refr_tmp<=idx.FF_test(end,nsub));
    idx.refr.FC.all(:,nsub) = refr_tmp(refr_tmp>=idx.FC_test(1,nsub) & refr_tmp<=idx.FC_test(end,nsub));
    
    %figure out the indices which denote the onset of a new block
    %(consider saving these in a cell array?)
    num_blocks = length(ics.blocksize{nsub});
    new_block_idx_tmp = new_block_dummy;
    new_block_idx_tmp(1:num_blocks) = cumsum(ics.blocksize{nsub})/2 + 1;
    idx.new_block_idx(:,nsub) = [1,new_block_idx_tmp];
    
    %combine all vEC trials together
    idx.vEC_all(:,nsub) = ec_tmp2;
    
    %segment refresher trials based on whether they were positive or
    %negative pertubations
    %NOTE:Update below to use tgt files!!
    
    if (ks==1 && strcmp(iraw.seq_id{nsub}, 'seq1')) ||...
            (ks==2 && strcmp(iraw.seq_id{nsub}, 'seq2'))
        %in this case, FF trials were surrounded by null trials, and FC
        %trials were surrounded by FC trials
        
        FC_refr_sign = sign(tgt_tmp(idx.refr.FC.all(:,nsub),6));
        FF_refr_sign = FC_refr_sign; %the correct thing to do would be to calculate the sign of the movement direction on the null trials
        
    elseif (ks==2 && strcmp(iraw.seq_id{nsub}, 'seq1')) ||...
            (ks==1 && strcmp(iraw.seq_id{nsub}, 'seq2'))
        %in this case, FF trials were surrounded by FC trials, and
        %FC trials were surrounded by null trials
        
        FF_refr_sign = sign(tgt_tmp(idx.refr.FF.all(:,nsub),6));
        FC_refr_sign = FF_refr_sign;
    end
    
    try
        idx.refr.FF.P(:, nsub) = idx.refr.FF.all(FF_refr_sign==1, nsub);    idx.refr.FF.N(:, nsub) = idx.refr.FF.all(FF_refr_sign==-1, nsub);
        idx.refr.FC.P(:, nsub) = idx.refr.FC.all(FC_refr_sign==1, nsub);    idx.refr.FC.N(:, nsub) = idx.refr.FC.all(FC_refr_sign==-1, nsub);
    catch
        keyboard;
    end
    
end

%get pre & post FF trials
idx.pre.FF.all = idx.pert.FF.all - 1;   idx.post.FF.all = idx.pert.FF.all + 1;
idx.pre.FF.P = idx.pert.FF.P - 1;   idx.post.FF.P = idx.pert.FF.P + 1;
idx.pre.FF.N = idx.pert.FF.N - 1;   idx.post.FF.N = idx.pert.FF.N + 1;

%combine training and test FC trials together
idx.FCp= [idx.FC_training; idx.FC_test];

%save each field's data
ics.idx = idx;
ics.dropped_FF_idx = dropped_FF_idx;
info.(cs) = ics;

%% save data for dat structure, and convert cell data to 3D matrix
%     dat.vxa_tmp.(cs) = vxa_tmp;
%     dat.vya_tmp.(cs) = vya_tmp;
%     dat.ypa_tmp.(cs) = ypa_tmp;
%     dat.xpa_tmp.(cs) = xpa_tmp;
%     dat.fra_tmp.(cs) = fra_tmp;
%     dat.xia_tmp.(cs) = xia_tmp;
%     dat.msa_tmp.(cs) = msa_tmp;
%     dat.fc_erra_tmp.(cs) = fc_erra_tmp;
%     dat.fKa_tmp.(cs) = fKa_tmp;
%     dat.fBa_tmp.(cs) = fBa_tmp;

%switch 3rd and 1st dimension to put data in trials x subject x sample format
dat.vxa_tmp.(cs) = permute(cell2mat(shiftdim(vxa_tmp,-2)),[3,2,1]);
dat.vya_tmp.(cs) = permute(cell2mat(shiftdim(vya_tmp,-2)),[3,2,1]);
dat.ypa_tmp.(cs) = permute(cell2mat(shiftdim(ypa_tmp,-2)),[3,2,1]);
dat.xpa_tmp.(cs) = permute(cell2mat(shiftdim(xpa_tmp,-2)),[3,2,1]);
dat.fra_tmp.(cs) = permute(cell2mat(shiftdim(fra_tmp,-2)),[3,2,1]);
dat.x_interp2a_tmp.(cs) = permute(cell2mat(shiftdim(x_interp2a_tmp,-2)),[3,2,1]);
dat.xp_bias_offseta_tmp.(cs) = permute(cell2mat(shiftdim(xp_bias_offseta_tmp,-2)),[3,2,1]);
dat.xp0a_tmp.(cs) = permute(cell2mat(shiftdim(xp0a_tmp,-2)),[3,2,1]);
dat.yp0a_tmp.(cs) = permute(cell2mat(shiftdim(yp0a_tmp,-2)),[3,2,1]);
dat.msa_tmp.(cs) = permute(cell2mat(shiftdim(msa_tmp,-2)),[3,2,1]);
dat.fera_tmp.(cs) = permute(cell2mat(shiftdim(fera_tmp,-2)),[3,2,1]);
dat.fKa_tmp.(cs) = permute(cell2mat(shiftdim(fKa_tmp,-2)),[3,2,1]);
dat.fBa_tmp.(cs) = permute(cell2mat(shiftdim(fBa_tmp,-2)),[3,2,1]);
dat.xr_ideala_tmp.(cs) = permute(cell2mat(shiftdim(xr_ideala_tmp,-2)),[3,2,1]);
dat.xs_ideala_tmp.(cs) = permute(cell2mat(shiftdim(xs_ideala_tmp,-2)),[3,2,1]);

dat.vymax.(cs) = vymax;
dat.tvymax.(cs) = tvymax;
dat.vtmax.(cs) = vtmax;
dat.tvtmax.(cs) = tvtmax;
dat.x0.(cs) = x0;
dat.y0.(cs) = y0;
dat.err_mm.(cs) = err_mm;
dat.err_deg.(cs) = err_deg;
dat.mov_dir.(cs) = mov_dir;
dat.Fcom.(cs) = Fcom;
dat.rt.(cs) = rt;
dat.mtime.(cs) = mtime;

info.num_samples = iraw.num_samples;
info.num_original_samples = iraw.num_samples;

%% determine and save indices based on the environment (FF/FC triplets in the LIE, and FF/FC triplets in the HIE)

% denv = extract_env_data_0827_2021a(dat, info);
denv = extract_env_data_0831_2021a(dat, info);
dat.denv = denv;

return


