%clear all;
close all;
home;


%jab did VZL and kam did ZLV


%% need to permutate the target files so that comparisons across subjects work
%in the future, data should be crunched with easy information about the exp
%sequence

%order them to be VZL

[total_trials, num_subjects] = size(dat.n);

exp_seq = {[1,2,3], [2,3,1]};

win = [-79:0];

num_samples = length(win);

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



