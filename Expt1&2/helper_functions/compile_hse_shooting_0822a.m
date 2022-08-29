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
original_samples = info.num_samples;
B = info.FFMAG; % 7.5 as of 8/27/2018

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
    sub_tgt_all_fwrd = sub_tgt_all(2:2:end,:);
    
    fwrd_tgt_all{k1} = sub_tgt_all_fwrd;
    tgt_all{k1} = sub_tgt_all;
    
    ec_tmp = find(sub_tgt_all_fwrd(:,7) > 0);
    idx.baseline_ec(k1,:) = ec_tmp(ec_tmp <= idx.baseline(k1,end));
    
    idx.vEC_training_ec(k1,:) = ec_tmp( ec_tmp<=idx.vEC.training(k1,end) & ec_tmp>=idx.vEC.training(k1,1) );
    idx.vEC_test_ec(k1,:) = ec_tmp( ec_tmp<=idx.vEC.test(k1,end) & ec_tmp>=idx.vEC.test(k1,1) );
    
    idx.zEC_training_ec(k1,:) = ec_tmp( ec_tmp<=idx.zEC.training(k1,end) & ec_tmp>=idx.zEC.training(k1,1) );
    idx.zEC_test_ec(k1,:) = ec_tmp( ec_tmp<=idx.zEC.test(k1,end) & ec_tmp>=idx.zEC.test(k1,1) );
    
    idx.null_training_ec(k1,:) = ec_tmp( ec_tmp<=idx.null.training(k1,end) & ec_tmp>=idx.null.training(k1,1) );
    idx.null_test_ec(k1,:) = ec_tmp( ec_tmp<=idx.null.test(k1,end) & ec_tmp>=idx.null.test(k1,1) );
    
    %find the triplets for each case
    triplets_tmp = find(sub_tgt_all_fwrd(:,2)~=0);
    
    vEC_triplets_tmp = triplets_tmp(triplets_tmp<=idx.vEC.test(k1,end) & triplets_tmp>=idx.vEC.test(k1,1) );
    zEC_triplets_tmp = triplets_tmp(triplets_tmp<=idx.zEC.test(k1,end) & triplets_tmp>=idx.zEC.test(k1,1) );
    null_triplets_tmp = triplets_tmp(triplets_tmp<=idx.null.test(k1,end) & triplets_tmp>=idx.null.test(k1,1) );
    
    idx.vEC_pre(k1,:) = vEC_triplets_tmp -1;     idx.vEC_post(k1,:) = vEC_triplets_tmp +1;
    idx.zEC_pre(k1,:) = zEC_triplets_tmp -1;     idx.zEC_post(k1,:) = zEC_triplets_tmp +1;
    idx.null_pre(k1,:) = null_triplets_tmp -1;     idx.null_post(k1,:) = null_triplets_tmp +1;
    
    
    %find the sign of the applied FF during the triplets
    sign_idx_vEC_tmp = vEC_triplets_tmp;
    sign_idx_zEC_tmp = zEC_triplets_tmp;
    sign_idx_null_tmp = null_triplets_tmp;
    
    FF_sign.vEC(k1,:) = sub_tgt_all_fwrd(sign_idx_vEC_tmp,2);
    FF_sign.zEC(k1,:) = sub_tgt_all_fwrd(sign_idx_zEC_tmp,2);
    FF_sign.null(k1,:) = sub_tgt_all_fwrd(sign_idx_null_tmp,2);
    
    %get all vEC trials
    vEC_tmp = find(sub_tgt_all_fwrd(:,7) > 0 & sub_tgt_all_fwrd(:,end) ~= 0);
    idx.vEC_all(k1,:) = vEC_tmp;
    
end

%% align data to the pillow

vel = cell(total_trials/2, 1);
ypos = cell(total_trials/2, 1);
fr = cell(total_trials/2, 1);

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
            
            if isempty(force_tmp)
                fr{cc}= [];
            else
                fr{cc}(:,k2) = force_tmp(win+pillow_idx, k2);
            end
            
        end
        
        cc=cc+1;
    end
end


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

frb = fr;


%% now calculate the AC on all trials

AC_all = nan(num_subjects, total_trials/2);

for k6=1:total_trials/2
    if ~isempty(frb{k6})
        for k7=1:num_subjects
            [b,~,~,~,s] = regress(frb{k6}(:,k7), [vel{k6}(:,k7)*0+1, vel{k6}(:,k7) * B]);
            AC_all(k7,k6) = b(2);
        end
    end
end

%% subtract baseline AC's

baseline_subtract_flag = 1;

for kq = 1:num_subjects
    
    AC.baseline(kq,:) = AC_all(kq, idx.baseline_ec(kq,:));
    
    if baseline_subtract_flag 
        %keyboard;
        AC_all(kq,:)  = AC_all(kq,:) - nanmean(AC.baseline(kq,:));
        baseline_subtract_flag = 0;
    end
    
end


%% flip the sign of the AC for the triplets
AC_unflipped_all = AC_all;
for k8=1:num_subjects
    
   AC_all(k8, idx.vEC_pre(k8,:)) = AC_all(k8, idx.vEC_pre(k8,:)) .* FF_sign.vEC(k8,:);
   AC_all(k8, idx.zEC_pre(k8,:)) = AC_all(k8, idx.zEC_pre(k8,:)) .* FF_sign.zEC(k8,:);
   AC_all(k8, idx.null_pre(k8,:)) = AC_all(k8, idx.null_pre(k8,:)) .* FF_sign.null(k8,:);
   
   AC_all(k8, idx.vEC_post(k8,:)) = AC_all(k8, idx.vEC_post(k8,:)) .* FF_sign.vEC(k8,:);
   AC_all(k8, idx.zEC_post(k8,:)) = AC_all(k8, idx.zEC_post(k8,:)) .* FF_sign.zEC(k8,:);
   AC_all(k8, idx.null_post(k8,:)) = AC_all(k8, idx.null_post(k8,:)) .* FF_sign.null(k8,:);
    
end

%% save the AC data in a struct

for kq=1:num_subjects
    AC.vEC.pre(kq,:) = AC_all(kq, idx.vEC_pre(kq,:));
    AC.zEC.pre(kq,:) = AC_all(kq, idx.zEC_pre(kq,:));
    AC.null.pre(kq,:) = AC_all(kq, idx.null_pre(kq,:));
    
    AC.vEC.post(kq,:) = AC_all(kq, idx.vEC_post(kq,:));
    AC.zEC.post(kq,:) = AC_all(kq, idx.zEC_post(kq,:));
    AC.null.post(kq,:) = AC_all(kq, idx.null_post(kq,:));
    
    AC.vEC.ar(kq,:) = AC.vEC.post(kq,:) - AC.vEC.pre(kq,:);
    AC.zEC.ar(kq,:) = AC.zEC.post(kq,:) - AC.zEC.pre(kq,:);
    AC.null.ar(kq,:) = AC.null.post(kq,:) - AC.null.pre(kq,:);

end
for kq=1:num_subjects
    AC_unflipped.vEC.pre(kq,:) = AC_unflipped_all(kq, idx.vEC_pre(kq,:));
    AC_unflipped.zEC.pre(kq,:) = AC_unflipped_all(kq, idx.zEC_pre(kq,:));
    AC_unflipped.null.pre(kq,:) = AC_unflipped_all(kq, idx.null_pre(kq,:));
    
    AC_unflipped.vEC.post(kq,:) = AC_unflipped_all(kq, idx.vEC_post(kq,:));
    AC_unflipped.zEC.post(kq,:) = AC_unflipped_all(kq, idx.zEC_post(kq,:));
    AC_unflipped.null.post(kq,:) = AC_unflipped_all(kq, idx.null_post(kq,:));
    
    AC_unflipped.vEC.ar(kq,:) = AC_unflipped.vEC.post(kq,:) - AC_unflipped.vEC.pre(kq,:);
    AC_unflipped.zEC.ar(kq,:) = AC_unflipped.zEC.post(kq,:) - AC_unflipped.zEC.pre(kq,:);
    AC_unflipped.null.ar(kq,:) = AC_unflipped.null.post(kq,:) - AC_unflipped.null.pre(kq,:);

end

%% plot the pre and post trials
figure;
LW=2;
x1 = [1:24];

for k=1:num_subjects,
    subplot(2,6,k); hold on;
    plot(AC.vEC.pre(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.vEC.pre(k,:))*ones(1,24), 'b');
    plot(AC.zEC.pre(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.zEC.pre(k,:))*ones(1,24), 'r');
    plot(AC.null.pre(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.null.pre(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,6,k+6); hold on;
    plot(AC_unflipped.vEC.pre(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.vEC.pre(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.zEC.pre(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.zEC.pre(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.null.pre(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.null.pre(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
end

suptitle('Pre AC trials');

figure;
for k=1:num_subjects,
    subplot(2,6,k); hold on;
    plot(AC.vEC.post(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.vEC.post(k,:))*ones(1,24), 'b');
    plot(AC.zEC.post(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.zEC.post(k,:))*ones(1,24), 'r');
    plot(AC.null.post(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.null.post(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,6,k+6); hold on;
    plot(AC_unflipped.vEC.post(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.vEC.post(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.zEC.post(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.zEC.post(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.null.post(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.null.post(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
end

suptitle('Post AC trials');

%% plot the adaptive response 


for k=1:num_subjects,
    subplot(2,6,k); hold on;
    plot(AC.vEC.ar(k,:), 'b', 'displayname', 'HSE (vEC)');  plot(x1,mean(AC.vEC.ar(k,:))*ones(1,24), 'b');
    plot(AC.zEC.ar(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC.zEC.ar(k,:))*ones(1,24), 'r');
    plot(AC.null.ar(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC.null.ar(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
    subplot(2,6,k+6); hold on;
    plot(AC_unflipped.vEC.ar(k,:), 'b', 'displayname', 'HSE (vEC)'); plot(x1,mean(AC_unflipped.vEC.ar(k,:))*ones(1,24), 'b');
    plot(AC_unflipped.zEC.ar(k,:), 'r', 'displayname', 'HSE (zEC)'); plot(x1,mean(AC_unflipped.zEC.ar(k,:))*ones(1,24), 'r');
    plot(AC_unflipped.null.ar(k,:), 'g', 'displayname', 'LSE (null)'); plot(x1,mean(AC_unflipped.null.ar(k,:))*ones(1,24), 'g');
    
    xlabel([info.sublist{k}, ': ', info.exp_seq{k}]);
    
end
for k=num_subjects+1,
    
    subplot(2,6,k); hold on;
    plot(x1, nanmean(AC.vEC.ar,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
    plot(x1, nanmean(AC.zEC.ar,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
    plot(x1, nanmean(AC.null.ar,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
    ylim([-2 2]);

    subplot(2,6,k+6); hold on;
    plot(x1, nanmean(AC_unflipped.vEC.ar,1), 'b', 'displayname', 'HSE (vEC)','linewidth',LW);
    plot(x1, nanmean(AC_unflipped.zEC.ar,1), 'r', 'displayname', 'HSE (zEC)','linewidth',LW);
    plot(x1, nanmean(AC_unflipped.null.ar,1), 'g', 'displayname', 'LSE (null)','linewidth',LW);
    ylim([-2 2]);
end
for k=1:12, subplot(2,6,k); ylim([-2 2]); end

suptitle('Adaptive response');

%% plot the force profiles, with the ideal, for each participant

%f_vEC_pre = get_sub_data(fbr, idx.vEC_pre);

%cfield = 'vEC';

for i=1:num_subjects,
    figure; 
    for k=1:24
        
        subplot(4,6,k);hold on;
        
        %plot the actual pre and post forces
        
        plot(frb{idx.vEC_pre(i,k)}(:,i), 'r');
        plot(frb{idx.vEC_post(i,k)}(:,i), 'b');
        
        %plot the ideal
        plot(vel{idx.vEC_pre(i,k)}(:,i) * B *-FF_sign.vEC(i,k) , 'k');
        
        
    end
    suptitle([info.sublist{i}, ': ', info.exp_seq{i}]);
    
end

%% plot the baseline force profiles


for i=1:num_subjects,
    figure; 
    for k=1:length(idx.baseline_ec(i,:))
        
        subplot(4,5,k);hold on;
        
        %plot the actual pre and post forces
        
        plot(frb{idx.baseline_ec(i,k)}(:,i), 'b');
        
        %plot the ideal
        plot(vel{idx.baseline_ec(i,k)}(:,i) * 7.5 , 'k');
        
        
    end
    suptitle([info.sublist{i}]);
    
end








%%



