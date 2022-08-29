clear all;
close all;
home;

cd('C:\Users\Laith Alhussein\Desktop\Research\HSE_exp\data\no_washout\shooting_data\mat_files\ALL_data');

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
hse_data=load(strcat(Cdr,'\',hse_filename));
dat=hse_data.dat;
info=hse_data.info;

[total_trials, num_subjects] = size(dat.n);
exp_seq = {[1,2,3], [2,3,1]};
win = [-79:0];
num_samples = length(win);
original_samples = info.num_samples;
B = info.FFMAG; % 7.5 as of 8/27/2018

flds = {'vEC', 'zEC', 'null'};

%% first permute the data to make indexing easier
    %we have VZL, ZLV, and LVZ
% for k1a = 1:num_subjects
%     
%     if strcmp(info.exp_seq{k1a}, 'ZLV')==1
%         dat = permute_data(dat, [3,1,2], k1a);
%     
%     elseif strcmp(info.exp_seq{k1a}, 'LVZ')==1
%         dat = permute_data(dat, [2,3,1], k1a);
%     end
%     
% end

%get the final tgt file



%% get the correct indices for each participant

%specifiy the blocks for each period, including both directions for now
% fam_block = 1;
% baseline_block = [2,3];
% training1_block = [4,5];
% test1_block = [6,7,8,9];
% training2_block = [10,11];
% test2_block = [12,13,14,15];
% training3_block = [16,17];
% test3_block = [18,19,20,21];

fam_idx = [1:50];
baseline_idx = [fam_idx(end)+ 1 : fam_idx(end) + 200];
training1_idx = [baseline_idx(end) + 1 : baseline_idx(end) + 202]; %the last few trials of the second training block is where the triplets start
test1_idx = [training1_idx(end) + 1 : training1_idx(end) +  516];
training2_idx = [test1_idx(end) + 1: test1_idx(end) + 202];
test2_idx = [training2_idx(end) + 1 : training2_idx(end) + 516];
training3_idx = [test2_idx(end) + 1: test2_idx(end) + 202];
test3_idx = [training3_idx(end) + 1: training3_idx(end) + 516];

%concentrate on forward trials...(?)
fam_idx = fam_idx(2:2:end)/2;
baseline_idx = baseline_idx(2:2:end)/2;
training1_idx = training1_idx(2:2:end)/2;
test1_idx = test1_idx(2:2:end)/2;
training2_idx = training2_idx(2:2:end)/2;
test2_idx = test2_idx(2:2:end)/2;
training3_idx = training3_idx(2:2:end)/2;
test3_idx = test3_idx(2:2:end)/2;

% fam_idx = [2:2:50];
% baseline_idx = [fam_idx(end)+ 2 :2: fam_idx(end) + 100];
% training1_idx = [baseline_idx(end) + 2 :2: baseline_idx(end) + 102]; %the last few trials of the second training block is where the triplets start
% test1_idx = [training1_idx(end) + 2 :2: training1_idx(end) +  258];
% training2_idx = [test1_idx(end) + 2:2: test1_idx(end) + 102];
% test2_idx = [training2_idx(end) + 2 :2: training2_idx(end) + 258];
% training3_idx = [test2_idx(end) + 2:2: test2_idx(end) + 102];
% test3_idx = [training3_idx(end) + 2:2: training3_idx(end) + 258];


tgt_all = cell(num_subjects,1);
fwrd_tgt_all = cell(num_subjects,1);

for k1=1:num_subjects
    %we have VZL, ZLV, and LVZ
    
    idx.baseline.all(k1,:) = baseline_idx;
    
    if strcmp(info.exp_seq{k1}, 'VZL')
       idx.training.vEC.all(k1,:) = training1_idx;
       idx.test.vEC.all(k1,:) = test1_idx;
        
       idx.training.zEC.all(k1,:) = training2_idx;
       idx.test.zEC.all(k1,:) = test2_idx;
       
       idx.training.null.all(k1,:) = training3_idx;
       idx.null.test.all(k1,:) = test3_idx;
       
    elseif strcmp(info.exp_seq{k1}, 'ZLV')
        idx.training.vEC.all(k1,:) = training3_idx;
        idx.test.vEC.all(k1,:) = test3_idx;
        
        idx.training.zEC.all(k1,:) = training1_idx;
        idx.test.zEC.all(k1,:) = test1_idx;
        
        idx.training.null.all(k1,:) = training2_idx;
        idx.null.test.all(k1,:) = test2_idx;
        
    elseif strcmp(info.exp_seq{k1}, 'LVZ')
        idx.training.vEC.all(k1,:) = training2_idx;
        idx.test.vEC.all(k1,:) = test2_idx;
        
        idx.training.zEC.all(k1,:) = training3_idx;
        idx.test.zEC.all(k1,:) = test3_idx;
        
        idx.training.null.all(k1,:) = training1_idx;
        idx.null.test.all(k1,:) = test1_idx;
    end
    
    sub_tgt_tmp = info.tgt_dat(k1,:);
    sub_tgt_all = vertcat(sub_tgt_tmp{:});
    sub_tgt_all_fwrd = sub_tgt_all(2:2:end,:);
    
    fwrd_tgt_all{k1} = sub_tgt_all_fwrd;
    tgt_all{k1} = sub_tgt_all_fwrd;
    
    ec_tmp = find(sub_tgt_all_fwrd(:,7) > 0);
    idx.baseline.ec(k1,:) = ec_tmp(ec_tmp <= idx.baseline.all(k1,end));
    
    idx.vEC.training.ec(k1,:) = ec_tmp( ec_tmp<=idx.training.vEC.all(k1,end) & ec_tmp>=idx.training.vEC.all(k1,1) );
    idx.vEC.test.ec(k1,:) = ec_tmp( ec_tmp<=idx.test.vEC.all(k1,end) & ec_tmp>=idx.test.vEC.all(k1,1) );
    
    idx.zEC.training.ec(k1,:) = ec_tmp( ec_tmp<=idx.training.zEC.all(k1,end) & ec_tmp>=idx.training.zEC.all(k1,1) );
    idx.zEC.test.ec(k1,:) = ec_tmp( ec_tmp<=idx.test.zEC.all(k1,end) & ec_tmp>=idx.test.zEC.all(k1,1) );
    
    idx.null.training.ec(k1,:) = ec_tmp( ec_tmp<=idx.training.null.all(k1,end) & ec_tmp>=idx.training.null.all(k1,1) );
    idx.null.test.ec(k1,:) = ec_tmp( ec_tmp<=idx.null.test.all(k1,end) & ec_tmp>=idx.null.test.all(k1,1) );
    
    %find the triplets for each case
    triplets_tmp = find(sub_tgt_all_fwrd(:,2)~=0);
    
    vEC_triplets_tmp = triplets_tmp(triplets_tmp<=idx.test.vEC.all(k1,end) & triplets_tmp>=idx.test.vEC.all(k1,1) );
    zEC_triplets_tmp = triplets_tmp(triplets_tmp<=idx.test.zEC.all(k1,end) & triplets_tmp>=idx.test.zEC.all(k1,1) );
    null_triplets_tmp = triplets_tmp(triplets_tmp<=idx.null.test.all(k1,end) & triplets_tmp>=idx.null.test.all(k1,1) );
    
    idx.vEC.pre(k1,:) = vEC_triplets_tmp -1;     idx.vEC.post(k1,:) = vEC_triplets_tmp +1;
    idx.zEC.pre(k1,:) = zEC_triplets_tmp -1;     idx.zEC.post(k1,:) = zEC_triplets_tmp +1;
    idx.null.pre(k1,:) = null_triplets_tmp -1;     idx.null.post(k1,:) = null_triplets_tmp +1;
    
    
    %find the sign of the applied FF during the triplets
    sign_idx_vEC_tmp = vEC_triplets_tmp;
    sign_idx_zEC_tmp = zEC_triplets_tmp;
    sign_idx_null_tmp = null_triplets_tmp;
    
    FF_sign.vEC(k1,:) = sub_tgt_all_fwrd(sign_idx_vEC_tmp,2);
    FF_sign.zEC(k1,:) = sub_tgt_all_fwrd(sign_idx_zEC_tmp,2);
    FF_sign.null(k1,:) = sub_tgt_all_fwrd(sign_idx_null_tmp,2);
    
    %get all vEC trials
    vEC_tmp = find(sub_tgt_all_fwrd(:,7) > 0 & sub_tgt_all_fwrd(:,end) ~= 0);
    idx.vEC.all(k1,:) = vEC_tmp;
    
end

%% align data to the pillow

vel = cell(total_trials/2, 1);
ypos = cell(total_trials/2, 1);
xpos = cell(total_trials/2, 1);
fr = cell(total_trials/2, 1);

err = nan(num_subjects, total_trials/2); %we should also calculate the error that the participant experienced

short_movements = 0;
junk_forces = 0;
for k2=1:num_subjects
    cc=1;
    for k3=2:2:total_trials
        
        pillow_idx = find( dat.pxr{k3}(:,k2) > 0.1, 1, 'first') - 1;
        
        force_tmp = dat.fyr{k3};
        
        %we can have cases in the early trials where there is no force data saved (these are just null trials anyway)
        if size(force_tmp,2)~=num_subjects, force_tmp(:,end:num_subjects) = NaN; end
        
        %sometimes the movement can be short, and we dont have enough data
        %before the pillow to subsample the movement
        if ~isempty(pillow_idx) & pillow_idx<length(win)
            vel{cc}(:,k2) = nan+win;
            ypos{cc}(:,k2) = nan+win;
            xpos{cc}(:,k2) = nan+win;
            if ~isempty(force_tmp)
                fr{cc}(:,k2) = nan+win;
            else
                fr{cc} = [];
            end
            short_movements = short_movements+1;
            
        else
            if isempty(pillow_idx), pillow_idx = abs(min(win)) + 1; end
            
            vel{cc}(:,k2) = dat.vxr{k3}(pillow_idx+win,k2);
            ypos{cc}(:,k2) = dat.pxr{k3}(pillow_idx+win,k2);
            xpos{cc}(:,k2) = dat.pyr{k3}(pillow_idx+win,k2);
            
            err(k2, cc) = dat.pyr{k3}(pillow_idx,k2) * 1000; %keep in mm

            if isempty(force_tmp)
                fr{cc}= [];
            else
                fr{cc}(:,k2) = force_tmp(win+pillow_idx, k2);
            end
            
        end
        
        cc=cc+1;
    end
end

%% now compare the data to when it is aligned to the midpoint of the movement

vel_mid = cell(total_trials/2, 1);
ypos_mid = cell(total_trials/2, 1);
xpos_mid = cell(total_trials/2, 1);
fr_mid = cell(total_trials/2, 1);

left_win = 60;
right_win = 20; %take 20 samples after the midpoint
total_win_mid = [-(left_win+right_win-1):0];

samples_to_ignore = 80; %ignore the first x samples when finding index at which the 5cm mark is passed

far_from_target = 0;

for k2=1:num_subjects
    cc=1;
    for k3=2:2:total_trials
        
        %if k3==26, keyboard; end
        
        pillow_idx = find( dat.pxr{k3}([samples_to_ignore+1:end],k2) > 0.1, 1, 'first') - 1 + samples_to_ignore;
        mid_idx =  find( dat.pxr{k3}([samples_to_ignore+1:end],k2) > 0.05, 1, 'first') - 1 + samples_to_ignore;
        %mid_idx_check = find( dat.pxr{k3}(:,k2) > 0.05 & dat.pxr{k3}(:,k2) < 0.06);
        
        force_tmp = dat.fyr{k3};
        
        %we can have cases in the early trials where there is no force data saved (these are just null trials anyway)
        if size(force_tmp,2)~=num_subjects, force_tmp(:,end:num_subjects) = NaN; end
        
        %sometimes the movement can be short, and we dont have enough data
        %before the pillow to subsample the movement
        if ~isempty(pillow_idx) & pillow_idx<length(win)
            vel_mid{cc}(:,k2) = nan+total_win_mid;
            ypos_mid{cc}(:,k2) = nan+total_win_mid;
            xpos_mid{cc}(:,k2) = nan+total_win_mid;
            if ~isempty(force_tmp)
                fr_mid{cc}(:,k2) = nan+total_win_mid;
            else
                fr_mid{cc} = [];
            end
            
        else
            if isempty(pillow_idx),
                pillow_idx = abs(min(total_win_mid)) + 1;
                mid_idx = pillow_idx;
                %mid_idx_check = mid_idx;
            end
            
%             if (mid_idx_check(end) - mid_idx_check(1)) * 0.005 > .5 %sometimes movements are already outside the starting target....
%                 far_from_target = far_from_target + 1;
%                 mid_idx = mid_idx_check(end) - 1;
%             
%             end
                
                try
                    vel_mid{cc}(:,k2) = dat.vxr{k3}([mid_idx - left_win: (mid_idx + right_win -1)],k2);
                catch
                    keyboard;
                end
                ypos_mid{cc}(:,k2) = dat.pxr{k3}([mid_idx - left_win: (mid_idx + right_win -1)],k2);
                xpos_mid{cc}(:,k2) = dat.pyr{k3}([mid_idx - left_win: (mid_idx + right_win -1)],k2);
                
                if isempty(force_tmp)
                    fr_mid{cc}= [];
                else
                    fr_mid{cc}(:,k2) = force_tmp([mid_idx - left_win: (mid_idx + right_win -1)], k2);
                end
                
                %if k3==26*2, keyboard; end
                
                %make sure to pad the data with nans after hitting the pillow
                if (mid_idx + right_win -1) > pillow_idx
                    num_pillow_samples = (mid_idx + right_win -1) - pillow_idx;
                    vel_mid{cc}(end-num_pillow_samples+1:end,k2) = NaN;
                    ypos_mid{cc}(end-num_pillow_samples+1,k2) = NaN;
                    xpos_mid{cc}(end-num_pillow_samples+1,k2) = NaN;
                    fr_mid{cc}(end-num_pillow_samples+1,k2) = NaN;
                    
                end
            
        end %pillow index is OK
        
        cc=cc+1;
    end
end

%keyboard;


%% subtract the baseline force profiles

%IMPORTANT: See Andrew's thesis to correctly estimate the adaptation on vEC
%trials. If we calculate it traditionally, then we are conflating learning
%with the effects of stiffness. But we can regress the AC's onto the
%direction of the vECs to estimate the stiffness effects, and then subtract
%those out. For now though, continue to subtract the baseline force
%profiles even for these trials, but make it optional in the future
% 
% vEC_bln_sub_flag = 1;
% 
% bln_fr = nan(num_subjects,length(win));
% bln_fr_tmp = get_sub_data(fr, idx.baseline_ec);
% 
% frb = cell(length(fr),1); %baseline subtracted force
% 
% for k4=1:num_subjects
%     
%     %bln_fr_tmp = get_sub_data(fr, idx.baseline_ec(k4,:));
%     %average across trials
%     bln_fr = squeeze(nanmean(bln_fr_tmp,2));
%     
%     for k5=1:length(vel)
%         if ~isempty(fr{k5})
%             
%             if ismember(k5,idx.vEC_all(k4,:))
%                 
%                 if vEC_bln_sub_flag, frb{k5}(:,k4) = fr{k5}(:,k4) - bln_fr(k4,:)'; end
%                 
%             else
%                 frb{k5}(:,k4) = fr{k5}(:,k4) - bln_fr(k4,:)';
%             end
%         end
%     end
% end

%frb = fr;


%% now calculate the AC on all trials, by regression and integration

AC_all = nan(num_subjects, total_trials/2);
AC_all_mid = AC_all;

ACI_all = AC_all;

ACI_all_mid = ACI_all;

for k6=1:total_trials/2
    if ~isempty(fr{k6})
        for k7=1:num_subjects
            [b,~,~,~,s] = regress(fr{k6}(:,k7), [vel{k6}(:,k7)*0+1, vel{k6}(:,k7) * B]);
            AC_all(k7,k6) = b(2);
            
            [b_mid,~,~,~,s_mid] = regress(fr_mid{k6}(:,k7), [vel_mid{k6}(:,k7)*0+1, vel_mid{k6}(:,k7) * B]);
            AC_all_mid(k7,k6) = b_mid(2);
            
            %ACI_all(k7,k6) = trapz(fr{k6}(:,k7)) / trapz(vel{k6}(:,k7) * B);
            ACI_all(k7,k6) = sum(fr{k6}(:,k7)) / sum( vel{k6}(:,k7) * B);
            ACI_all_mid(k7,k6) = sum(fr_mid{k6}(:,k7)) / sum( vel_mid{k6}(:,k7) * B);
            
        end
    end
end


%% subtract baseline AC's

baseline_subtract_flag = 1;

for kq = 1:num_subjects
    
    AC.baseline(kq,:) = AC_all(kq, idx.baseline.ec(kq,:));
    AC_mid.baseline(kq,:) = AC_all_mid(kq, idx.baseline.ec(kq,:));
    
    ACI.baseline(kq,:) = ACI_all(kq, idx.baseline.ec(kq,:));
    ACI_mid.baseline(kq,:) = ACI_all_mid(kq, idx.baseline.ec(kq,:));
        
    if baseline_subtract_flag 
        %keyboard;
        AC_all(kq,:)  = AC_all(kq,:) - nanmean(AC.baseline(kq,:));
        AC_all_mid(kq,:)  = AC_all_mid(kq,:) - nanmean(AC_mid.baseline(kq,:));
        
        ACI_all(kq,:) = ACI_all(kq,:) - nanmean(ACI.baseline(kq,:));
        ACI_all_mid(kq,:) = ACI_all_mid(kq,:) - nanmean(ACI_mid.baseline(kq,:));
        
        baseline_subtract_flag = 0;
    end
    
end


%% flip the sign of the AC for the triplets
AC_unflipped_all = AC_all;
AC_unflipped_all_mid = AC_all_mid;
for k8=1:num_subjects
    
   AC_all(k8, idx.vEC.pre(k8,:)) = AC_all(k8, idx.vEC.pre(k8,:)) .* -FF_sign.vEC(k8,:);
   AC_all(k8, idx.zEC.pre(k8,:)) = AC_all(k8, idx.zEC.pre(k8,:)) .* -FF_sign.zEC(k8,:);
   AC_all(k8, idx.null.pre(k8,:)) = AC_all(k8, idx.null.pre(k8,:)) .* -FF_sign.null(k8,:);
   
   AC_all(k8, idx.vEC.post(k8,:)) = AC_all(k8, idx.vEC.post(k8,:)) .* -FF_sign.vEC(k8,:);
   AC_all(k8, idx.zEC.post(k8,:)) = AC_all(k8, idx.zEC.post(k8,:)) .* -FF_sign.zEC(k8,:);
   AC_all(k8, idx.null.post(k8,:)) = AC_all(k8, idx.null.post(k8,:)) .* -FF_sign.null(k8,:);
   
   
   %%%repeat for midpoint case
   AC_all_mid(k8, idx.vEC.pre(k8,:)) = AC_all_mid(k8, idx.vEC.pre(k8,:)) .* -FF_sign.vEC(k8,:);
   AC_all_mid(k8, idx.zEC.pre(k8,:)) = AC_all_mid(k8, idx.zEC.pre(k8,:)) .* -FF_sign.zEC(k8,:);
   AC_all_mid(k8, idx.null.pre(k8,:)) = AC_all_mid(k8, idx.null.pre(k8,:)) .* -FF_sign.null(k8,:);
   
   AC_all_mid(k8, idx.vEC.post(k8,:)) = AC_all_mid(k8, idx.vEC.post(k8,:)) .* -FF_sign.vEC(k8,:);
   AC_all_mid(k8, idx.zEC.post(k8,:)) = AC_all_mid(k8, idx.zEC.post(k8,:)) .* -FF_sign.zEC(k8,:);
   AC_all_mid(k8, idx.null.post(k8,:)) = AC_all_mid(k8, idx.null.post(k8,:)) .* -FF_sign.null(k8,:);
   
   
   ACI_all(k8, idx.vEC.pre(k8,:)) = ACI_all(k8, idx.vEC.pre(k8,:)) .* -FF_sign.vEC(k8,:);
   ACI_all(k8, idx.zEC.pre(k8,:)) = ACI_all(k8, idx.zEC.pre(k8,:)) .* -FF_sign.zEC(k8,:);
   ACI_all(k8, idx.null.pre(k8,:)) = ACI_all(k8, idx.null.pre(k8,:)) .* -FF_sign.null(k8,:);
   
   ACI_all(k8, idx.vEC.post(k8,:)) = ACI_all(k8, idx.vEC.post(k8,:)) .* -FF_sign.vEC(k8,:);
   ACI_all(k8, idx.zEC.post(k8,:)) = ACI_all(k8, idx.zEC.post(k8,:)) .* -FF_sign.zEC(k8,:);
   ACI_all(k8, idx.null.post(k8,:)) = ACI_all(k8, idx.null.post(k8,:)) .* -FF_sign.null(k8,:);
   
   
   %ACI_all_mid
   ACI_all_mid(k8, idx.vEC.pre(k8,:)) = ACI_all_mid(k8, idx.vEC.pre(k8,:)) .* -FF_sign.vEC(k8,:);
   ACI_all_mid(k8, idx.zEC.pre(k8,:)) = ACI_all_mid(k8, idx.zEC.pre(k8,:)) .* -FF_sign.zEC(k8,:);
   ACI_all_mid(k8, idx.null.pre(k8,:)) = ACI_all_mid(k8, idx.null.pre(k8,:)) .* -FF_sign.null(k8,:);
   
   ACI_all_mid(k8, idx.vEC.post(k8,:)) = ACI_all_mid(k8, idx.vEC.post(k8,:)) .* -FF_sign.vEC(k8,:);
   ACI_all_mid(k8, idx.zEC.post(k8,:)) = ACI_all_mid(k8, idx.zEC.post(k8,:)) .* -FF_sign.zEC(k8,:);
   ACI_all_mid(k8, idx.null.post(k8,:)) = ACI_all_mid(k8, idx.null.post(k8,:)) .* -FF_sign.null(k8,:);
   
   
end

%% save the AC data in a struct

for kq=1:num_subjects
    AC.pre.vEC(kq,:) = AC_all(kq, idx.vEC.pre(kq,:));
    AC.pre.zEC(kq,:) = AC_all(kq, idx.zEC.pre(kq,:));
    AC.pre.null(kq,:) = AC_all(kq, idx.null.pre(kq,:));
    
    AC.post.vEC(kq,:) = AC_all(kq, idx.vEC.post(kq,:));
    AC.post.zEC(kq,:) = AC_all(kq, idx.zEC.post(kq,:));
    AC.post.null(kq,:) = AC_all(kq, idx.null.post(kq,:));
    
    AC.ar.vEC(kq,:) = AC.post.vEC(kq,:) - AC.pre.vEC(kq,:);
    AC.ar.zEC(kq,:) = AC.post.zEC(kq,:) - AC.pre.zEC(kq,:);
    AC.ar.null(kq,:) = AC.post.null(kq,:) - AC.pre.null(kq,:);
    
    %%%repeat for mid case
    AC_mid.pre.vEC(kq,:) = AC_all_mid(kq, idx.vEC.pre(kq,:));
    AC_mid.pre.zEC(kq,:) = AC_all_mid(kq, idx.zEC.pre(kq,:));
    AC_mid.pre.null(kq,:) = AC_all_mid(kq, idx.null.pre(kq,:));
    
    AC_mid.post.vEC(kq,:) = AC_all_mid(kq, idx.vEC.post(kq,:));
    AC_mid.post.zEC(kq,:) = AC_all_mid(kq, idx.zEC.post(kq,:));
    AC_mid.post.null(kq,:) = AC_all_mid(kq, idx.null.post(kq,:));
    
    AC_mid.ar.vEC(kq,:) = AC_mid.post.vEC(kq,:) - AC.pre.vEC(kq,:);
    AC_mid.ar.zEC(kq,:) = AC_mid.post.zEC(kq,:) - AC.pre.zEC(kq,:);
    AC_mid.ar.null(kq,:) = AC_mid.post.null(kq,:) - AC.pre.null(kq,:);

    %%%unflipped
    AC_unflipped.pre.vEC(kq,:) = AC_unflipped_all(kq, idx.vEC.pre(kq,:));
    AC_unflipped.pre.zEC(kq,:) = AC_unflipped_all(kq, idx.zEC.pre(kq,:));
    AC_unflipped.pre.null(kq,:) = AC_unflipped_all(kq, idx.null.pre(kq,:));
    
    AC_unflipped.post.vEC(kq,:) = AC_unflipped_all(kq, idx.vEC.post(kq,:));
    AC_unflipped.post.zEC(kq,:) = AC_unflipped_all(kq, idx.zEC.post(kq,:));
    AC_unflipped.post.null(kq,:) = AC_unflipped_all(kq, idx.null.post(kq,:));
    
    AC_unflipped.ar.vEC(kq,:) = AC_unflipped.post.vEC(kq,:) - AC_unflipped.pre.vEC(kq,:);
    AC_unflipped.ar.zEC(kq,:) = AC_unflipped.post.zEC(kq,:) - AC_unflipped.pre.zEC(kq,:);
    AC_unflipped.ar.null(kq,:) = AC_unflipped.post.null(kq,:) - AC_unflipped.pre.null(kq,:);
    
    
    %%%unflipped for mid case
    AC_mid_unflipped.pre.vEC(kq,:) = AC_unflipped_all_mid(kq, idx.vEC.pre(kq,:));
    AC_mid_unflipped.pre.zEC(kq,:) = AC_unflipped_all_mid(kq, idx.zEC.pre(kq,:));
    AC_mid_unflipped.pre.null(kq,:) = AC_unflipped_all_mid(kq, idx.null.pre(kq,:));
    
    AC_mid_unflipped.post.vEC(kq,:) = AC_unflipped_all_mid(kq, idx.vEC.post(kq,:));
    AC_mid_unflipped.post.zEC(kq,:) = AC_unflipped_all_mid(kq, idx.zEC.post(kq,:));
    AC_mid_unflipped.post.null(kq,:) = AC_unflipped_all_mid(kq, idx.null.post(kq,:));
    
    AC_mid_unflipped.ar.vEC(kq,:) = AC_mid_unflipped.post.vEC(kq,:) - AC_mid_unflipped.pre.vEC(kq,:);
    AC_mid_unflipped.ar.zEC(kq,:) = AC_mid_unflipped.post.zEC(kq,:) - AC_mid_unflipped.pre.zEC(kq,:);
    AC_mid_unflipped.ar.null(kq,:) = AC_mid_unflipped.post.null(kq,:) - AC_mid_unflipped.pre.null(kq,:);
    
    
    %%%repeat for the integration
    ACI.pre.vEC(kq,:) = ACI_all(kq, idx.vEC.pre(kq,:));
    ACI.pre.zEC(kq,:) = ACI_all(kq, idx.zEC.pre(kq,:));
    ACI.pre.null(kq,:) = ACI_all(kq, idx.null.pre(kq,:));
    
    ACI.post.vEC(kq,:) = ACI_all(kq, idx.vEC.post(kq,:));
    ACI.post.zEC(kq,:) = ACI_all(kq, idx.zEC.post(kq,:));
    ACI.post.null(kq,:) = ACI_all(kq, idx.null.post(kq,:));
    
    ACI.ar.vEC(kq,:) = ACI.post.vEC(kq,:) - ACI.pre.vEC(kq,:);
    ACI.ar.zEC(kq,:) = ACI.post.zEC(kq,:) - ACI.pre.zEC(kq,:);
    ACI.ar.null(kq,:) = ACI.post.null(kq,:) - ACI.pre.null(kq,:);
    
        %%%integration for the mid case
    ACI.pre.vEC(kq,:) = ACI_all(kq, idx.vEC.pre(kq,:));
    ACI.pre.zEC(kq,:) = ACI_all(kq, idx.zEC.pre(kq,:));
    ACI.pre.null(kq,:) = ACI_all(kq, idx.null.pre(kq,:));
    
    ACI.post.vEC(kq,:) = ACI_all(kq, idx.vEC.post(kq,:));
    ACI.post.zEC(kq,:) = ACI_all(kq, idx.zEC.post(kq,:));
    ACI.post.null(kq,:) = ACI_all(kq, idx.null.post(kq,:));
    
    ACI.ar.vEC(kq,:) = ACI.post.vEC(kq,:) - ACI.pre.vEC(kq,:);
    ACI.ar.zEC(kq,:) = ACI.post.zEC(kq,:) - ACI.pre.zEC(kq,:);
    ACI.ar.null(kq,:) = ACI.post.null(kq,:) - ACI.pre.null(kq,:);
    
    
    %get the AC for the vEC periods
    ACI.training.vEC(kq,:) = ACI_all(kq, idx.training.vEC.all(kq,:));
    ACI.test.vEC(kq,:) = ACI_all(kq, idx.test.vEC.all(kq,:));
    
    
    
end

%% plot the pre and post trials
figure;
LW=2;
x1 = [1:24];

for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(AC.pre.vEC(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.pre.vEC(k,:))*ones(1,24), 'b');
    plot(AC.pre.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.pre.zEC(k,:))*ones(1,24), 'r');
    plot(AC.pre.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.pre.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(AC_unflipped.pre.vEC(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.pre.vEC(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.pre.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.pre.zEC(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.pre.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.pre.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
end

suptitle('Pre AC trials');

figure;
for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(AC.post.vEC(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.post.vEC(k,:))*ones(1,24), 'b');
    plot(AC.post.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.post.zEC(k,:))*ones(1,24), 'r');
    plot(AC.post.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.post.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(AC_unflipped.post.vEC(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.post.vEC(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.post.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.post.zEC(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.post.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.post.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
end

suptitle('Post AC trials');

%% plot the adaptive response 

figure;
for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(AC.ar.vEC(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.ar.vEC(k,:))*ones(1,24), 'b');
    plot(AC.ar.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.ar.zEC(k,:))*ones(1,24), 'r');
    plot(AC.ar.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.ar.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(AC_unflipped.ar.vEC(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.ar.vEC(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.ar.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.ar.zEC(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.ar.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.ar.null(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
end

%figure;
for k=1:num_subjects,
    
    subplot(2,num_subjects,k); hold on;
    plot(x1, nanmean(AC.ar.vEC,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
    plot(x1, nanmean(AC.ar.zEC,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
    plot(x1, nanmean(AC.ar.null,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
    ylim([-2 2]);

    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(x1, nanmean(AC_unflipped.ar.vEC,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
    plot(x1, nanmean(AC_unflipped.ar.zEC,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
    plot(x1, nanmean(AC_unflipped.ar.null,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
    ylim([-2 2]);
end
for k=1:12, subplot(2,num_subjects,k); ylim([-2 2]); end

suptitle('Adaptive response');

%% plot the force profiles, with the ideal, for each participant

%f_vEC_pre = get_sub_data(fbr, idx.vEC_pre);

%cfield = 'vEC';

for i=1:num_subjects,
    figure; 
    for k=1:24
        
        subplot(4,num_subjects,k);hold on;
        
        %plot the actual pre and post forces
        
        plot(fr{idx.null.pre(i,k)}(:,i), 'r');
        plot(fr{idx.null.post(i,k)}(:,i), 'b');
        
        %plot the ideal
        plot(vel{idx.null.pre(i,k)}(:,i) * B *-FF_sign.vEC(i,k) , 'k');
        
        %plot the AR
        plot(fr{idx.null.post(i,k)}(:,i) - fr{idx.null.pre(i,k)}(:,i), 'color', [0.5, 0, 0.5]);
        
        grid on;
        
        title(['ACI is: ', num2str(ACI.ar.null(i,k),2)] );
        
    end
    suptitle(['Pre & Post: ', info.sublist{i}, ': ', info.exp_seq{i}]);
    
end

%% plot the force profiles, with the ideal, for each participant

%f_vEC_pre = get_sub_data(fbr, idx.vEC_pre);                               
                                                                                                      
%% plot the baseline force profiles


% for i=1:num_subjects,
%     figure; 
%     for k=1:length(idx.baseline.ec(i,:))
%         
%         subplot(4,num_subjects,k); hold on;
%         
%         %plot the actual pre and post forces
%         
%         plot(fr_mid{idx.baseline.ec(i,k)}(:,i), 'b');
%         
%         %plot the ideal
%         plot(vel_mid{idx.baseline.ec(i,k)}(:,i) * B , 'k');
%         
%         %plot the position
%         plot(ypos_mid{idx.baseline.ec(i,k)}(:,i) * 100);
%         
%         grid on; grid minor;
%         
%         
%     end
%     suptitle(['Baseline, ',info.sublist{i}]);
%     
% end

%% calculate the average FP for each case

f_sub.pre.vEC = get_sub_data(fr, idx.vEC.pre);
f_sub.pre.zEC = get_sub_data(fr, idx.zEC.pre);
f_sub.pre.null = get_sub_data(fr, idx.null.pre);

f_sub.post.vEC = get_sub_data(fr, idx.vEC.post);
f_sub.post.zEC = get_sub_data(fr, idx.zEC.post);
f_sub.post.null = get_sub_data(fr, idx.null.post);

f_sub.ar.vEC = f_sub.post.vEC - f_sub.pre.vEC;
f_sub.ar.zEC = f_sub.post.zEC - f_sub.pre.zEC;
f_sub.ar.null = f_sub.post.null - f_sub.pre.null;

vel_sub.pre.vEC = get_sub_data(vel, idx.vEC.pre);
vel_sub.pre.zEC = get_sub_data(vel, idx.zEC.pre);
vel_sub.pre.null = get_sub_data(vel, idx.null.pre);

vel_sub.post.vEC = get_sub_data(vel, idx.vEC.post);
vel_sub.post.zEC = get_sub_data(vel, idx.zEC.post);
vel_sub.post.null = get_sub_data(vel, idx.null.post);

ypos_sub.pre.vEC = get_sub_data(ypos, idx.vEC.pre);
ypos_sub.pre.zEC = get_sub_data(ypos, idx.zEC.pre);
ypos_sub.pre.null = get_sub_data(ypos, idx.null.pre);

ypos_sub.post.vEC = get_sub_data(ypos, idx.vEC.post);
ypos_sub.post.zEC = get_sub_data(ypos, idx.zEC.post);
ypos_sub.post.null = get_sub_data(ypos, idx.null.post);


%%%repeat for when the movement is aligned by the middle
fmid_sub.pre.vEC = get_sub_data(fr_mid, idx.vEC.pre);
fmid_sub.pre.zEC = get_sub_data(fr_mid, idx.zEC.pre);
fmid_sub.pre.null = get_sub_data(fr_mid, idx.null.pre);

fmid_sub.post.vEC = get_sub_data(fr_mid, idx.vEC.post);
fmid_sub.post.zEC = get_sub_data(fr_mid, idx.zEC.post);
fmid_sub.post.null = get_sub_data(fr_mid, idx.null.post);

fmid_sub.ar.vEC = fmid_sub.post.vEC - fmid_sub.pre.vEC;
fmid_sub.ar.zEC = fmid_sub.post.zEC - fmid_sub.pre.zEC;
fmid_sub.ar.null = fmid_sub.post.null - fmid_sub.pre.null;

velmid_sub.pre.vEC = get_sub_data(vel_mid, idx.vEC.pre);
velmid_sub.pre.zEC = get_sub_data(vel_mid, idx.zEC.pre);
velmid_sub.pre.null = get_sub_data(vel_mid, idx.null.pre);

velmid_sub.post.vEC = get_sub_data(vel_mid, idx.vEC.post);
velmid_sub.post.zEC = get_sub_data(vel_mid, idx.zEC.post);
velmid_sub.post.null = get_sub_data(vel_mid, idx.null.post);

yposmid_sub.pre.vEC = get_sub_data(ypos_mid, idx.vEC.pre);
yposmid_sub.pre.zEC = get_sub_data(ypos_mid, idx.zEC.pre);
yposmid_sub.pre.null = get_sub_data(ypos_mid, idx.null.pre);

yposmid_sub.post.vEC = get_sub_data(ypos_mid, idx.vEC.post);
yposmid_sub.post.zEC = get_sub_data(ypos_mid, idx.zEC.post);
yposmid_sub.post.null = get_sub_data(ypos_mid, idx.null.post);

%% flip the sign of the FF for the adaptive response

f2_sub.ar.vEC = f_sub.ar.vEC * NaN;
f2_sub.ar.zEC = f_sub.ar.zEC * NaN;
f2_sub.ar.null = f_sub.ar.null * NaN;

f2mid_sub.ar.vEC = fmid_sub.ar.vEC * NaN;
f2mid_sub.ar.zEC = fmid_sub.ar.zEC * NaN;
f2mid_sub.ar.null = fmid_sub.ar.null * NaN;

for q1=1:num_subjects
    
    %if q1==4, keyboard; end
    
    for q2 = 1:24
        f2_sub.ar.vEC(q1,q2,:) = f_sub.ar.vEC(q1,q2,:) * -FF_sign.vEC(q1,q2);
        f2_sub.ar.zEC(q1,q2,:) = f_sub.ar.zEC(q1,q2,:) * -FF_sign.zEC(q1,q2);
        f2_sub.ar.null(q1,q2,:) = f_sub.ar.null(q1,q2,:) * -FF_sign.null(q1,q2);
        
        f2mid_sub.ar.vEC(q1,q2,:) = fmid_sub.ar.vEC(q1,q2,:) * -FF_sign.vEC(q1,q2);
        f2mid_sub.ar.zEC(q1,q2,:) = fmid_sub.ar.zEC(q1,q2,:) * -FF_sign.zEC(q1,q2);
        f2mid_sub.ar.null(q1,q2,:) = fmid_sub.ar.null(q1,q2,:) * -FF_sign.null(q1,q2);
        
    end
    
    %the below method cant handle NaNs in matrix mutltiplications
%     f2_sub.ar.vEC(q1,:,:) = permute( shiftdim( squeeze(f_sub.ar.vEC(q1,:,:))' * diag(-FF_sign.vEC(q1,:)), -1), [1, 3, 2]);
%     f2_sub.ar.zEC(q1,:,:) = permute( shiftdim( squeeze(f_sub.ar.zEC(q1,:,:))' * diag(-FF_sign.zEC(q1,:)), -1), [1, 3, 2]);
%     f2_sub.ar.null(q1,:,:) = permute( shiftdim( squeeze(f_sub.ar.null(q1,:,:))' * diag(-FF_sign.null(q1,:)), -1), [1, 3, 2]);
    
end
%% plot the AC for the entire vEC period for each participant

%for each subject, make a 3 x 3 subplot (experiment on y-axis, pre, post,
%and AR on x-axis) for averaged force profiles
%each panel should have the force profiles with different alignments, and
%the ideal for each case

purple = [0.5, 0, 0.5];
LW2 = 1.5;

% SE_flag = 0;
% 
% for kk=1:num_subjects
%     
%     %force_vEC_tmp = get_sub_data(fr, 
%     
%     figure;
%     
%     %first is vEC
%     subplot(3,3,1); hold on;
%     plot( squeeze( nanmean(f_sub.pre.vEC(kk,:,:), 2) ), 'r', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.pre.vEC(kk,:,:), 2) ), 'r--', 'linewidth', LW2);
%     
%     
%     plot( squeeze( nanmean(vel_sub.pre.vEC(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.pre.vEC(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.pre.vEC(kk,:,:), 2), nanstd(f_sub.pre.vEC(kk,:,:), 0, 2),...
%             [1:size(f_sub.pre.vEC,3)], size(f_sub.pre.vEC,2),'r');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.pre.vEC(kk,:,:), 2), nanstd(fmid_sub.pre.vEC(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.pre.vEC,3)], size(fmid_sub.pre.vEC,2),'r');
%     end
%     
%     %next is post
%     subplot(3,3,2); hold on;
%     plot( squeeze( nanmean(f_sub.post.vEC(kk,:,:), 2) ), 'b', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.post.vEC(kk,:,:), 2) ), 'b--', 'linewidth', LW2);
%     
%     plot( squeeze( nanmean(vel_sub.post.vEC(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.post.vEC(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.post.vEC(kk,:,:), 2), nanstd(f_sub.post.vEC(kk,:,:), 0, 2),...
%             [1:size(f_sub.post.vEC,3)], size(f_sub.post.vEC,2),'b');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.post.vEC(kk,:,:), 2), nanstd(fmid_sub.post.vEC(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.post.vEC,3)], size(fmid_sub.post.vEC,2),'b');
%     end
%     
%     
%     %adaptive response
%     subplot(3,3,3); hold on;
%     plot( squeeze( nanmean(f2_sub.ar.vEC(kk,:,:), 2) ), 'color', purple, 'linewidth', LW2);
%     plot( squeeze( nanmean(f2mid_sub.ar.vEC(kk,:,:), 2) ), 'color', purple, 'Linestyle', '--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f2_sub.ar.vEC(kk,:,:), 2), nanstd(f2_sub.ar.vEC(kk,:,:), 0, 2),...
%             [1:size(f2_sub.ar.vEC,3)], size(f2_sub.ar.vEC,2),purple);
%         standard_error_shading_07_16_2015(nanmean(f2mid_sub.ar.vEC(kk,:,:), 2), nanstd(f2mid_sub.ar.vEC(kk,:,:), 0, 2),...
%             [1:size(f2mid_sub.ar.vEC,3)], size(f2mid_sub.ar.vEC,2),purple);
%     end
%     
%    %next do the zEC case
%    %pre
%     subplot(3,3,4); hold on;
%     plot( squeeze( nanmean(f_sub.pre.zEC(kk,:,:), 2) ), 'r', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.pre.zEC(kk,:,:), 2) ), 'r--', 'linewidth', LW2);
%     
%     plot( squeeze( nanmean(vel_sub.pre.zEC(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.pre.zEC(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.pre.zEC(kk,:,:), 2), nanstd(f_sub.pre.zEC(kk,:,:), 0, 2),...
%             [1:size(f_sub.pre.zEC,3)], size(f_sub.pre.zEC,2),'r');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.pre.zEC(kk,:,:), 2), nanstd(fmid_sub.pre.zEC(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.pre.zEC,3)], size(fmid_sub.pre.zEC,2),'r');
%     end
%     
%     %next is post
%     subplot(3,3,5); hold on;
%     plot( squeeze( nanmean(f_sub.post.zEC(kk,:,:), 2) ), 'b', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.post.zEC(kk,:,:), 2) ), 'b--', 'linewidth', LW2);
%     
%     plot( squeeze( nanmean(vel_sub.post.zEC(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.post.zEC(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.post.zEC(kk,:,:), 2), nanstd(f_sub.post.zEC(kk,:,:), 0, 2),...
%             [1:size(f_sub.post.zEC,3)], size(f_sub.post.zEC,2),'b');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.post.zEC(kk,:,:), 2), nanstd(fmid_sub.post.zEC(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.post.zEC,3)], size(fmid_sub.post.zEC,2),'b');
%     end
%     
%     %adaptive response
%     subplot(3,3,6); hold on;
%     plot( squeeze( nanmean(f2_sub.ar.zEC(kk,:,:), 2) ), 'color', purple, 'linewidth', LW2);
%     plot( squeeze( nanmean(f2mid_sub.ar.zEC(kk,:,:), 2) ), 'color', purple, 'Linestyle', '--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f2_sub.ar.zEC(kk,:,:), 2), nanstd(f2_sub.ar.zEC(kk,:,:), 0, 2),...
%             [1:size(f2_sub.ar.zEC,3)], size(f2_sub.ar.zEC,2),purple);
%         standard_error_shading_07_16_2015(nanmean(f2mid_sub.ar.zEC(kk,:,:), 2), nanstd(f2mid_sub.ar.zEC(kk,:,:), 0, 2),...
%             [1:size(f2mid_sub.ar.zEC,3)], size(f2mid_sub.ar.zEC,2),purple);
%     end
%     
%     %%%next do the null case
%     
%     %first is vEC
%     subplot(3,3,7); hold on;
%     plot( squeeze( nanmean(f_sub.pre.null(kk,:,:), 2) ), 'r', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.pre.null(kk,:,:), 2) ), 'r--', 'linewidth', LW2);
%     
%     
%     plot( squeeze( nanmean(vel_sub.pre.null(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.pre.null(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.pre.null(kk,:,:), 2), nanstd(f_sub.pre.null(kk,:,:), 0, 2),...
%             [1:size(f_sub.pre.null,3)], size(f_sub.pre.null,2),'r');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.pre.null(kk,:,:), 2), nanstd(fmid_sub.pre.null(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.pre.null,3)], size(fmid_sub.pre.null,2),'r');
%     end
%     
%     %next is post
%     subplot(3,3,8); hold on;
%     plot( squeeze( nanmean(f_sub.post.null(kk,:,:), 2) ), 'b', 'linewidth', LW2);
%     plot( squeeze( nanmean(fmid_sub.post.null(kk,:,:), 2) ), 'b--', 'linewidth', LW2);
%     
%     plot( squeeze( nanmean(vel_sub.post.null(kk,:,:), 2) ) * B, 'k', 'linewidth', LW2);
%     plot( squeeze( nanmean(velmid_sub.post.null(kk,:,:), 2) ) * B, 'k--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f_sub.post.null(kk,:,:), 2), nanstd(f_sub.post.null(kk,:,:), 0, 2),...
%             [1:size(f_sub.post.null,3)], size(f_sub.post.null,2),'b');
%         standard_error_shading_07_16_2015(nanmean(fmid_sub.post.null(kk,:,:), 2), nanstd(fmid_sub.post.null(kk,:,:), 0, 2),...
%             [1:size(fmid_sub.post.null,3)], size(fmid_sub.post.null,2),'b');
%     end
%     
%     %adaptive response
%     subplot(3,3,9); hold on;
%     plot( squeeze( nanmean(f2_sub.ar.null(kk,:,:), 2) ), 'color', purple, 'linewidth', LW2);
%     plot( squeeze( nanmean(f2mid_sub.ar.null(kk,:,:), 2) ), 'color', purple, 'Linestyle', '--', 'linewidth', LW2);
%     
%     if SE_flag
%         standard_error_shading_07_16_2015(nanmean(f2_sub.ar.null(kk,:,:), 2), nanstd(f2_sub.ar.null(kk,:,:), 0, 2),...
%             [1:size(f2_sub.ar.null,3)], size(f2_sub.ar.null,2),purple);
%         standard_error_shading_07_16_2015(nanmean(f2mid_sub.ar.null(kk,:,:), 2), nanstd(f2mid_sub.ar.null(kk,:,:), 0, 2),...
%             [1:size(f2mid_sub.ar.null,3)], size(f2mid_sub.ar.null,2),purple);
%     end
%     
%     suptitle([info.sublist{kk}]);
% end

%% to figure out which alignment to go with, we should calculate a variance metric over a window, and compare them

w2 = [25:65]; %window over which we will calculate the variance for each case

%we will do this for pre and post trials separately
for i=1:length(flds)
    for kk=1:num_subjects
        cfld = flds{i};
        fs.pre.(cfld)(kk,1) =  std(nanmean(f_sub.pre.(cfld)(kk,:,w2), 2)); %signal
        fs.pre.(cfld)(kk,2) =  mean(nanstd(f_sub.pre.(cfld)(kk,:,w2),0,2)); %noise
        
        fs.post.(cfld)(kk,1) = std(nanmean( f_sub.post.(cfld)(kk,:,w2),2));
        fs.post.(cfld)(kk,2) = mean(nanstd(f_sub.post.(cfld)(kk,:,w2),0,2));
        
        fsmid.pre.(cfld)(kk,1) = std(nanmean( fmid_sub.pre.(cfld)(kk,:,w2),2));
        fsmid.pre.(cfld)(kk,2) = mean(nanstd(fmid_sub.pre.(cfld)(kk,:,w2),0,2));
        
        fsmid.post.(cfld)(kk,1) = std(nanmean( fmid_sub.post.(cfld)(kk,:,w2),2));
        fsmid.post.(cfld)(kk,2) = mean(nanstd(fmid_sub.post.(cfld)(kk,:,w2),0,2));
    end
end

%3 x 2 subplot ( signal, noise, ratio vs pre and post)
nn=[1:num_subjects];

figure;

%flds2 = {'pre', 'post'};
%i+(k1-1)*6

%subplot(3,2,kq + (kk-1)*2); hold on;

subplot(321); hold on; title('Signal, pre');

%signal first
plot(fs.pre.vEC(:,1), 'b'); plot(fsmid.pre.vEC(:,1), 'b--');
plot(fs.pre.zEC(:,1), 'r'); plot(fsmid.pre.zEC(:,1), 'r--');
plot(fs.pre.null(:,1), 'g'); plot(fsmid.pre.null(:,1), 'g--');

subplot(322); hold on; title('Signal, post');

plot(fs.post.vEC(:,1), 'b'); plot(fsmid.post.vEC(:,1), 'b--');
plot(fs.post.zEC(:,1), 'r'); plot(fsmid.post.zEC(:,1), 'r--');
plot(fs.post.null(:,1), 'g'); plot(fsmid.post.null(:,1), 'g--');

%noise
subplot(323); hold on; title('Noise, pre');

plot(fs.pre.vEC(:,2), 'b'); plot(fsmid.pre.vEC(:,2), 'b--');
plot(fs.pre.zEC(:,2), 'r'); plot(fsmid.pre.zEC(:,2), 'r--');
plot(fs.pre.null(:,2), 'g'); plot(fsmid.pre.null(:,2), 'g--');

subplot(324); hold on; title('Noise, post');

plot(fs.post.vEC(:,2), 'b'); plot(fsmid.post.vEC(:,2), 'b--');
plot(fs.post.zEC(:,2), 'r'); plot(fsmid.post.zEC(:,2), 'r--');
plot(fs.post.null(:,2), 'g'); plot(fsmid.post.null(:,2), 'g--');

%SNR
subplot(325); hold on; title('snr, pre');
plot(fs.pre.vEC(:,1) ./ fs.pre.vEC(:,2), 'b'); plot(fsmid.pre.vEC(:,1) ./ fsmid.pre.vEC(:,2), 'b--');
plot(fs.pre.zEC(:,1) ./ fs.pre.zEC(:,2), 'r'); plot(fsmid.pre.zEC(:,1) ./ fsmid.pre.zEC(:,2), 'r--');
plot(fs.pre.null(:,1) ./ fs.pre.null(:,2), 'g'); plot(fsmid.pre.null(:,1) ./ fsmid.pre.null(:,2), 'g--');

subplot(326); hold on; title('snr, post');
plot(fs.post.vEC(:,1) ./ fs.post.vEC(:,2), 'b'); plot(fsmid.post.vEC(:,1) ./ fsmid.post.vEC(:,2), 'b--');
plot(fs.post.zEC(:,1) ./ fs.post.zEC(:,2), 'r'); plot(fsmid.post.zEC(:,1) ./ fsmid.post.zEC(:,2), 'r--');
plot(fs.post.null(:,1) ./ fs.post.null(:,2), 'g'); plot(fsmid.post.null(:,1) ./ fsmid.post.null(:,2), 'g--');


%% we want to see if we have an effect immediately after training, so focus on first 5 or so trials, and look at the force profiles

%make a 3 x 3 plot, of pre, post and ar, and experiment on y axis

%num_early_trials = 4;

% for kk=1:num_subjects
%     figure;
%     
%     for jj=1:num_early_trials
%         
%         %vEC
%        subplot(331); hold on;
%        plot(squeeze(f_sub.pre.vEC(kk,jj,:)), 'r');
%        
%        subplot(332); hold on;
%        plot(squeeze(f_sub.post.vEC(kk,jj,:)), 'b');
%        
%        subplot(333); hold on;
%        plot(squeeze(f_sub.post.vEC(kk,jj,:)), 'color', purple);
%        
%        %zEC
%        subplot(334); hold on;
%        plot(squeeze(f_sub.pre.zEC(kk,jj,:)), 'r');
%        
%        subplot(335); hold on;
%        plot(squeeze(f_sub.post.zEC(kk,jj,:)), 'b');
%        
%        subplot(336); hold on;
%        plot(squeeze(f_sub.post.zEC(kk,jj,:)), 'color', purple);
%        
%        %null
%        subplot(337); hold on;
%        plot(squeeze(f_sub.pre.null(kk,jj,:)), 'r');
%        
%        subplot(338); hold on;
%        plot(squeeze(f_sub.post.null(kk,jj,:)), 'b');
%        
%        subplot(339); hold on;
%        plot(squeeze(f_sub.post.null(kk,jj,:)), 'color', purple);
%         
%     end
%     
%     subplot(331); hold on;
%     plot(squeeze(nanmean(f_sub.pre.vEC(kk,[1:num_early_trials],:))), 'r', 'linewidth', 2);
%     
%     subplot(332); hold on;
%     plot(squeeze(nanmean(f_sub.post.vEC(kk,[1:num_early_trials],:))), 'b', 'linewidth', 2);
%     
%     subplot(333); hold on;
%     plot(squeeze(nanmean(f_sub.ar.vEC(kk,[1:num_early_trials],:))), 'color', purple, 'linewidth', 2);
%     
%     
%     subplot(334); hold on;
%     plot(squeeze(nanmean(f_sub.pre.zEC(kk,[1:num_early_trials],:))), 'r', 'linewidth', 2);
%     
%     subplot(335); hold on;
%     plot(squeeze(nanmean(f_sub.post.zEC(kk,[1:num_early_trials],:))), 'b', 'linewidth', 2);
%     
%     subplot(336); hold on;
%     plot(squeeze(nanmean(f_sub.ar.zEC(kk,[1:num_early_trials],:))), 'color', purple, 'linewidth', 2);
%     
%     
%     subplot(337); hold on;
%     plot(squeeze(nanmean(f_sub.pre.null(kk,[1:num_early_trials],:))), 'r', 'linewidth', 2);
%     
%     subplot(338); hold on;
%     plot(squeeze(nanmean(f_sub.post.null(kk,[1:num_early_trials],:))), 'b', 'linewidth', 2);
%     
%     subplot(339); hold on;
%     plot(squeeze(nanmean(f_sub.ar.null(kk,[1:num_early_trials],:))), 'color', purple, 'linewidth', 2);
%     
%     suptitle([ info.sublist{kk}, ': ', info.exp_seq{kk}]);
%     
% 
% end


%% instead of including pre and post, just make 1 plot with all subjects and plot the ar
%make sure to add the AC...

num_early_trials = 5;

figure;

for p1 = 1:3
    for kk=1:num_subjects
        
        cfld = flds{p1};
        
        %first do vEC
        subplot(3,num_subjects, kk + (p1-1) * num_subjects); hold on;
        for jj=1:num_early_trials
            plot(squeeze(smooth( f2_sub.ar.(cfld)(kk,jj,:), 7) ) , 'color', purple);
        end
        
        %show the mean AR and ideal (get ideal from pre trials)
        plot( squeeze( nanmean(f2_sub.ar.(cfld)(kk, [1:num_early_trials], :), 2) ), 'color', purple, 'linewidth', 2);
        plot( squeeze( nanmean(vel_sub.pre.(cfld)(kk, [1:num_early_trials], :) * B, 2) ), 'k', 'linewidth', 2);
        
        ylim([-8 8]); 
        
        if p1==1, title([info.sublist{kk}, ': ', info.exp_seq{kk}]); end
        
        if kk==1
            if p1==1
                ylabel(cfld);
            elseif p1==2
                ylabel(cfld);
            elseif p1==3
                ylabel(cfld);
            end
        end
        
    end
    
end

suptitle('Early Adaptive Response');

%keyboard;

%% look at learning curves during test period for vEC

for k=1:num_subjects,
    subplot(1,num_subjects,k); hold on;
    plot(ACI.ar.vEC(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,nanmean(ACI.ar.vEC(k,:))*ones(1,24), 'b');
    plot(ACI.ar.zEC(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,nanmean(ACI.ar.zEC(k,:))*ones(1,24), 'r');
    plot(ACI.ar.null(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,nanmean(ACI.ar.null(k,:))*ones(1,24), 'g');
    
%     plot(x1, nanmean(ACI.ar.vEC,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
%     plot(x1, nanmean(ACI.ar.zEC,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
%     plot(x1, nanmean(ACI.ar.null,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
    %ylim([-10 10]);
end
%for k=1:12, subplot(1,num_subjects,k); ylim([-2 2]); end

suptitle('Adaptive response by integration');


%% look at learning curves during test period for vEC


figure; hold on;
plot(nanmean(ACI.ar.vEC,1), 'b', 'displayname', 'HSE (vEC)');  plot(x1,nanmean(ACI.ar.vEC(:))*ones(1,24), 'b');
plot(nanmean(ACI.ar.zEC,1), 'r', 'displayname', 'HSE (zEC)'); plot(x1,nanmean(ACI.ar.zEC(:))*ones(1,24), 'r');
plot(nanmean(ACI.ar.null,1), 'g', 'displayname', 'LSE (null)'); plot(x1,nanmean(ACI.ar.null(:))*ones(1,24), 'g');

plot(cumsum(nanmean(ACI.ar.vEC,1))./[1:24], 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
plot(cumsum(nanmean(ACI.ar.zEC,1))./[1:24], 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
plot(cumsum(nanmean(ACI.ar.null,1))./[1:24], 'g', 'displayname', 'LSE (null)','linewidth',LW);

%     plot(x1, nanmean(ACI.ar.vEC,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
%     plot(x1, nanmean(ACI.ar.zEC,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
%     plot(x1, nanmean(ACI.ar.null,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
%ylim([-10 10]);
%for k=1:12, subplot(1,num_subjects,k); ylim([-2 2]); end

suptitle('Adaptive response by integration - mean across subjects');


%% first look to see how well the imposed errors correlate with the experienced errors

%get the imposed error sequence in mm (separate based on training only or training + test)
ierr_tmp = tgt_all{1}(idx.vEC.all(1,:), end); %this is in degrees
ierr.all  = tand(ierr_tmp) *100; %in mm

%repeat ierr since we will combine participant data
ierr_rep = repmat(ierr.all, [num_subjects,1]);

%concatinate the participant data
err_vEC_all = nan(num_subjects * length(ierr_tmp), 1);

c1 = 1;
for k=1:num_subjects
    for j=1:length(ierr_tmp)
        err_vEC_all(c1) = err(k, idx.vEC.all(k,j));
        c1=c1+1;
    end
end

%plot the experienced (y) vs imposed (x)

figure; hold on;
subplot(121); hold on;
plot(ierr_rep, err_vEC_all, '.');
ylabel('experienced error');
xlabel('imposed error');
%calculate the correlation
[corr_e, p_e] = corr(ierr_rep, err_vEC_all, 'rows', 'complete');
title(['r = ', num2str(corr_e), ' and p = ', num2str(p_e)]);
axis([-20, 20, -20, 20]);

subplot(122);
plot(ierr_rep, err_vEC_all - ierr_rep,'.');
ylabel('experienced - imposed error');
xlabel('imposed error');
title(['\sigma = ', num2str(nanstd(err_vEC_all - ierr_rep),3), ' mm']);
xlim([-20 20])


%keyboard;
%% look at the AC, and do the modeling for the vEC training period

burn_trials = 10; %a certain number of trials need to be removed to exclude the effects of the EC bias
ACI_training_vEC_cut = ACI.training.vEC(:, [burn_trials+1:end]);


%subtract each participant's mean response to mitigate EC bias
EC_bias_flag = 1;
if EC_bias_flag
   ACI_training_vEC_cut = bsxfun(@(x,y) x-y, ACI_training_vEC_cut, nanmean(ACI_training_vEC_cut,2));
    EC_bias_flag = 0;
end



ierr.training = ierr_tmp(1: length(idx.training.vEC.all(1,:)));
ierr_training_cut = ierr.training([burn_trials+1:end]);

%concatinate
ierr_training_rep = repmat(ierr_training_cut, [num_subjects,1]);
% ACI_training_vEC_all = reshape(ACI.training.vEC', [size(ACI.training.vEC,2) * size(ACI.training.vEC,1), 1]);
ACI_training_vEC_all = reshape(ACI_training_vEC_cut', [size(ACI_training_vEC_cut,2) * size(ACI_training_vEC_cut,1), 1]);

%to estimate the effect of stiffness, regress the adpatation onto the
%imposed error encountered on the same trial

%do the regression
X1 = [ierr_training_rep*0+1, ierr_training_rep];
[b2,~,~,~,s2] = regress(ACI_training_vEC_all, X1);


ACI_avg_training_vEC = nanmean(ACI_training_vEC_cut, 1);

figure; hold on;
plot(ierr_training_rep, ACI_training_vEC_all, '.');
%axis([-10, 10, -2.5, 2.5]);

figure; hold on;
plot(ierr_training_cut, ACI_avg_training_vEC, '.');

% for k=1:num_subjects
%     plot(ierr.training, ACI.training.vEC(k,:), '.');
% 
% end

corr(ierr_training_rep, ACI_training_vEC_all, 'rows', 'complete')

corr(ierr_training_cut, ACI_avg_training_vEC')

axis([-10, 10, -1, 1]);

%% bootstrap for model fitting
% n_boot = 5000;
% A = nan(n_boot,1);
% B = nan(n_boot,1);
% C = nan(n_boot,1);
% 
% for j=1:n_boot
%     csel = randsample(num_subjects, size(ACI_training_vEC_cut,2), 'true');
%     yboot_sample = diag(ACI_training_vEC_cut(csel,:));
%     xboot_sample = ierr_training_cut(csel);
%     qq = nlinfit(xboot_sample, yboot_sample, @stiffness_model, [0.1,0.1,0.1]);
%     A(j) = qq(1);
%     B(j) = qq(2);
%     C(j) = qq(3);
% end
% 
% figure; hist(A);
% figure; hist(B);
% figure; hist(C);

qq = nlinfit(ierr_training_cut, ACI_avg_training_vEC, @stiffness_model, [0.8,0.07,0.07]);


%%% remove effects of stiffness
ACI_avg_training_ns = ACI_avg_training_vEC - ( (0*ACI_avg_training_vEC+b2(1)) + (ierr_training_cut*b2(2))' ); %offset is first in b2

%fit the model
qq2 = nlinfit(ierr_training_cut, ACI_avg_training_ns, @stiffness_model2, [0.8,0.07,0.07]);

keyboard;





