clear all;
close all;
home;

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
hse_data=load(strcat(Cdr,'\',hse_filename));
dat=hse_data.dat;
info=hse_data.info;

[total_trials, num_subjects] = size(dat.n);
exp_seq = {[1,2,3], [2,3,1]};
win = [-79:0];
num_samples = length(win);

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
fam_block = 1;
baseline_block = [2,3];
training1_block = [4,5];
test1_block = [6,7,8,9];
training2_block = [10,11];
test2_block = [12,13,14,15];
training3_block = [16,17];
test3_block = [18,19,20,21];

fam_idx = [1:50];
baseline_idx = [fam_idx(end)+ 1 : fam_idx(end)+ 1 + 200];
training1_idx = [baseline_idx(end) + 1 : baseline_idx(end) + 1 + 208];
test1_idx = [training1_idx(end) + 1 : training1_idx(end) + 1 +  510];
training2_idx = [test1_idx(end) + 1: test1_idx(end) + 1 + 208];
test2_idx = [training2_idx(end) + 1 : training2_idx(end) + 1 + 510];
training3_idx = [test2_idx(end) + 1: test2_idx(end) + 1 + 208];
test3_idx = [training3_idx(end) + 1: training3_idx(end) + 1 + 510];

tgt_all = cell(num_subjects,1);

for k1=1:num_subjects
    %we have VZL, ZLV, and LVZ
    
    idx.baseline(k1,:) = baseline_idx;
    
    if strcmp(info.exp_seq{k1}, 'VZL')
       idx.vEC.training(k1,:) = training1_idx;
       idx.vEC.test(k1,:) = test1_idx;
        
       idx.zEC.training(k1,:) = training2_idx;
       idx.zEC.test(k1,:) = test2_idx;
       
       idx.null.training(k1,:) = training3_idx;
       idx.null.test(k1,:) = test3_idx;
       
    elseif strcmp(info.exp_seq{k1}, 'ZLV')
        idx.vEC.training(k1,:) = training3_idx;
        idx.vEC.test(k1,:) = test3_idx;
        
        idx.zEC.training(k1,:) = training1_idx;
        idx.zEC.test(k1,:) = test1_idx;
        
        idx.null.training(k1,:) = training2_idx;
        idx.null.test(k1,:) = test2_idx;
        
    elseif strcmp(info.exp_seq{k1}, 'LVZ')
        idx.vEC.training(k1,:) = training2_idx;
        idx.vEC.test(k1,:) = test2_idx;
        
        idx.zEC.training(k1,:) = training3_idx;
        idx.zEC.test(k1,:) = test3_idx;
        
        idx.null.training(k1,:) = training1_idx;
        idx.null.test(k1,:) = test1_idx;
    end
    
    sub_tgt_tmp = info.tgt_dat(k1,:);
    sub_tgt_all = vertcat(sub_tgt_tmp{:});
    tgt_all{k1} = sub_tgt_all;
    
    ec_tmp = find(sub_tgt_all(:,7) > 0);
    idx.baseline_ec = ec_tmp(ec_tmp <= idx.baseline(k1,end));
    
    idx.vEC_training_ec = ec_tmp( ec_tmp<=idx.vEC.training(k1,end) & ec_tmp>=idx.vEC.training(k1,1) );
    idx.vEC_test_ec = ec_tmp( ec_tmp<=idx.vEC.test(k1,end) & ec_tmp>=idx.vEC.test(k1,1) );
    
    idx.zEC_training_ec = ec_tmp( ec_tmp<=idx.zEC.training(k1,end) & ec_tmp>=idx.zEC.training(k1,1) );
    idx.zEC_test_ec = ec_tmp( ec_tmp<=idx.zEC.test(k1,end) & ec_tmp>=idx.zEC.test(k1,1) );
    
    idx.null_training_ec = ec_tmp( ec_tmp<=idx.null.training(k1,end) & ec_tmp>=idx.null.training(k1,1) );
    idx.null_test_ec = ec_tmp( ec_tmp<=idx.null.test(k1,end) & ec_tmp>=idx.null.test(k1,1) );
    
end

%% align data to the pillow

vel = cell(total_trials/2, 1);
ypos = cell(total_trials/2, 1);
fr = cell(total_trials/2, 1);

short_movements = 0;
for k2=1:num_subjects
    cc=1;
    for k3=2:2:total_trials
        
        pillow_idx = find( dat.pxr{k3}(:,k2) > 0.1, 1, 'first') - 1;
        
        %sometimes the movement can be short, and we dont have enough data
        %before the pillow to subsample the movement
        if ~isempty(pillow_idx) & pillow_idx<length(win)
            vel{cc}(:,k2) = nan+win;
            ypos{cc}(:,k2) = nan+win;
            fr{cc}(:,k2) = nan+win;
            
        else
            if isempty(pillow_idx), pillow_idx = abs(min(win)) + 1; end
            
            vel{cc}(:,k2) = dat.vxr{k3}(pillow_idx+win,k2);
            ypos{cc}(:,k2) = dat.pxr{k3}(pillow_idx+win,k2);
            
            force_tmp = dat.fyr{k3};
            if isempty(force_tmp)
                fr{cc}= [];
            else
                fr{cc}(:,k2) = force_tmp(win+pillow_idx, k2);
            end
            
        end
        
        cc=cc+1;
    end
end







keyboard;

%% need to permutate the target files so that comparisons across subjects work
%in the future, data should be crunched with easy information about the exp
%sequence

%order them to be VZL


vEC1_all = nan(num_subjects, 24,num_samples);
vEC2_all = nan(num_subjects, 24,num_samples);
vel_vEC1_all = nan(num_subjects, 24,num_samples);
vel_vEC2_all = nan(num_subjects, 24,num_samples);
ypos_vEC1_all = nan(num_subjects, 24,num_samples);
ypos_vEC2_all = nan(num_subjects, 24,num_samples);

ypos_vEC1_all_tmp = nan(num_subjects, 24,451);
ypos_vEC2_all_tmp = nan(num_subjects, 24,451);

FF_sign_vEC =nan(num_subjects, 24);


for k=1:num_subjects
    tgt_tmp = info.tgt_dat(1,:);
    sub_tgt_all = vertcat(tgt_tmp{:});
    
    %find the vEC period
    vEC_idx_first = find(sub_tgt_all(:,7)==2 & sub_tgt_all(:,9)~=0,1,'first');
    vEC_idx_last = find(sub_tgt_all(:,7)==2 & sub_tgt_all(:,9)~=0 ,1,'last') + 6; %offset 3 trials for the last triplet
    
    %focus on EC trials during this period (this will only come from the test period)
    EC_idx_all= find(sub_tgt_all(:,7)==1);
    vEC_force_idx = EC_idx_all(EC_idx_all >= vEC_idx_first & EC_idx_all <= vEC_idx_last);
    
    %now calculate the adaptive response in terms of the force
    
    first_EC = vEC_force_idx(1:2:end);
    last_EC = vEC_force_idx(2:2:end);
    
    
    for i=1:24%there are 24 triplets
        
        %find point of passing target for each case
        %get position
        ypos_vEC1_all_tmp(k,i,:) = dat.pxr{first_EC(i)}(:,k);
        ypos_vEC2_all_tmp(k,i,:) = dat.pxr{last_EC(i)}(:,k);
        
        target_idx1 = find(ypos_vEC1_all_tmp(k,i,:) > 0.1,1,'first')-1;
        target_idx2 = find(ypos_vEC2_all_tmp(k,i,:) > 0.1,1,'first')-1;
        
%         if isempty(target_idx1), target_idx1 = 0; end
%         if isempty(target_idx2), target_idx2 = 0; end
        
        
        if ~isempty(target_idx1)
            vEC1_all(k,i,:) = dat.fyr{first_EC(i)}(target_idx1+win,k);
            vel_vEC1_all(k,i,:) = dat.vxr{first_EC(i)}(target_idx1+win,k);
            ypos_vEC1_all(k,i,:) = dat.pxr{first_EC(i)}(target_idx1+win,k);
        end
        
        if ~isempty(target_idx2)
            vEC2_all(k,i,:) = dat.fyr{last_EC(i)}(target_idx2+win,k);
            vel_vEC2_all(k,i,:) = dat.vxr{last_EC(i)}(target_idx2+win,k);
            ypos_vEC2_all(k,i,:) = dat.pxr{last_EC(i)}(target_idx2+win,k);
        end
        
        
        FF_sign_vEC(k,i) = sub_tgt_all(first_EC(i)+2,2);
        
    end
    
    
%     for kk=1:24
%        %get the sign of the FF
%         FF_sign_vEC(k,kk) = sub_tgt_all(first_EC(kk)+2,2); 
%     end
    
    
end






[AC_vEC1_all, rsq] = calculate_AC(vEC1_all, 7.5*FF_sign_vEC(1,:), vel_vEC1_all);
[AC_vEC2_all, rsq] = calculate_AC(vEC2_all, 7.5*FF_sign_vEC(1,:), vel_vEC2_all);


AR_all = AC_vEC2_all - AC_vEC1_all;






for i=1:num_subjects,
    figure;
    for k=1:24
        
        subplot(4,6,k);hold on;
        plot(squeeze(vEC1_all(i,k,:)), 'b'); 
         plot(squeeze(vEC2_all(i,k,:)), 'r');

        %plot(squeeze(vel_vEC2_all(i,k,:)) * 7.5, 'r');
        
        %plot(squeeze(ypos_vEC2_all(i,k,:)), 'g'); 
        
    end
end



