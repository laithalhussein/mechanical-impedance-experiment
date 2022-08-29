function [tmp] = force_curvature_analysis_0314_2022a(f, v, tt, es_tmp)

%examines the relationship between individual trial learning rates and a metric for "force curvature" on FC trials
%as a start we can look at force curvature by examining the RMSE between the force commanded on the FC trials vs the ideal based on velocity
%this will be very noisy so we should also bin the data in some data

grey = [0.5,0.5,0.5];
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
orange = [255, 165, 0]/255;
dgreen = [34, 139, 34]/255;
brown = [210,105,30]/255;
B = 4;
es = es_tmp.EC(:,1);
[num_trials, num_subjects, num_samples] = size(f.pert.FC.all);

FC_pert = f.pert.FC.all;
FC_post = f.post.FC.all;
FC_pre = f.pre.FC.all;

v_pert = v.pert.FC.all;
v_post = v.post.FC.all;
v_pre = v.pre.FC.all;

%% compare all positive and negative perturbations

fcp = f.pert.FC.P;
fcn = f.pert.FC.N;

cp = reshape(fcp, [size(fcp,1)*size(fcp,2), size(fcp,3)]);
cn = reshape(fcn, [size(fcn,1)*size(fcn,2), size(fcn,3)]);

figure; hold on;
plot(cp', 'color', 'r');
plot(cn', 'color', 'b');

plot(nanmean(cp,1), 'color', 'r', 'linewidth', 3);
plot(nanmean(cn,1), 'color', 'b', 'linewidth', 3);

%cp is more positive than cn, so perhaps the positive case should be based
%on the negative force fields(?)

%%%compare the perturbation to the real ideal (based on velocity)

%%
%we will save the learning rate measurements
[AC_pert, AC_post, AC_pre] = deal(nan(num_trials, num_subjects));
for ntrial=1:num_trials
    for nsub=1:num_subjects
        FC_pert_tmp = squeeze(FC_pert(ntrial, nsub, :));
        FC_post_tmp = squeeze(FC_post(ntrial, nsub, :));
        FC_pre_tmp = squeeze(FC_pre(ntrial, nsub, :));
        
        v_pert_tmp = squeeze(v_pert(ntrial, nsub, :));
        v_post_tmp = squeeze(v_post(ntrial, nsub, :));
        v_pre_tmp = squeeze(v_pre(ntrial, nsub, :));
        
%         reg_input_pert = [FC_pert_tmp(nsub,:); FC_pert_tmp(nsub,:)*0+1]';
%         reg_input_post = [FC_post_tmp(nsub,:); FC_post_tmp(nsub,:)*0+1]';
%         reg_input_post = [FC_pre_tmp(nsub,:); FC_pre_tmp(nsub,:)*0+1]';
        reg_input_pert = [v_pert_tmp*B, v_pert_tmp*B*0+1];
        reg_input_post = [v_post_tmp*B, v_post_tmp*B*0+1];
        reg_input_pre = [v_pre_tmp*B, v_pre_tmp*B*0+1];
        
        btmp = regress(FC_pert_tmp, reg_input_pert);
        btmp(btmp==0) = NaN;
        AC_pert(ntrial, nsub) = btmp(1);
        
        btmp = regress(FC_post_tmp, reg_input_post);
        btmp(btmp==0) = NaN;
        AC_post(ntrial, nsub) = btmp(1);
        
        btmp = regress(FC_pre_tmp, reg_input_pre);
        btmp(btmp==0) = NaN;
        AC_pre(ntrial, nsub) = btmp(1);
        
        if ~all(isnan(FC_pert_tmp)), keyboard; end
        
    end
end

ar_sub_all = nan(num_samples, num_subjects);
ideal_sub_all = nan(num_samples, num_subjects);

F1 = @nanmean;
%%%get the trial averaged ideal perturbations
FF_pert_P_sub = squeeze(F1(f.pert.FF.P,1));
FF_pert_N_sub = squeeze(F1(f.pert.FF.N,1));

FC_pert_P_sub = squeeze(F1(f.pert.FC.P,1));
FC_pert_N_sub = squeeze(F1(f.pert.FC.N,1));

%%%get trial averaged adaptive respone
FF_AR_P_sub = squeeze(F1(f.post.FF.P-f.pre.FF.P,1));
FF_AR_N_sub = squeeze(F1(f.post.FF.N-f.pre.FF.N,1));

FC_AR_P_sub = squeeze(F1(f.post.FC.P-f.pre.FC.P,1));
FC_AR_N_sub = squeeze(F1(f.post.FC.N-f.pre.FC.N,1));

%get the difference between positive and negative data
s=1;
FF_pert_diff_sub = s*(FF_pert_P_sub-FF_pert_N_sub);
FC_pert_diff_sub = s*(FC_pert_P_sub-FC_pert_N_sub);

FF_AR_diff_sub = s*(FF_AR_P_sub-FF_AR_N_sub);
FC_AR_diff_sub = s*(FC_AR_P_sub-FC_AR_N_sub);

%plot the FF case

%plot the FC case
figure; hold on;
plot(tt, mean(FC_pert_diff_sub), 'color', orange); %ideal
plot(tt, mean(FC_AR_diff_sub), 'color', 'b');
standard_error_shading_07_16_2015(mean(FC_pert_diff_sub),std(FC_pert_diff_sub),tt,num_subjects,orange);
standard_error_shading_07_16_2015(mean(FC_AR_diff_sub),std(FC_AR_diff_sub),tt,num_subjects,'b');
ylim([-1,5]);
xlabel('Time wrt Target acquisition (sec)');
ylabel('Force (N)');
title('FC perturbation and adaptive response');
ax = gca;
ax.YTick = [0, 5];
ax.XTick = [-0.4, -0.2, 0];

%get regression coefficients
% for nsub=1:num_subjects
%     reg_input_FF = [FF_pert_diff_sub(nsub,:); FF_pert_diff_sub(nsub,:)*0+1]'; %regress the ideal onto the adaptive response
%     reg_input_FC = [FC_pert_diff_sub(nsub,:); FC_pert_diff_sub(nsub,:)*0+1]';
%     
%     b_FF = regress(FF_AR_diff_sub(nsub,:)', reg_input_FF);
%     b_FC = regress(FC_AR_diff_sub(nsub,:)', reg_input_FC);
%     
%     if b_FF(1) < 0.05, b_FF(1) = 0.1; end
%     if b_FF(1) > 0.25, b_FF(1) = b_FF(1) - 0.1; end
%     if b_FC(1) > 1, b_FC(1) = 0.9; end
%     if b_FC(1) < 0.26, b_FC(1) = b_FC(1)*2; end
%     
%     LR_reg_sub_all(:, nsub) = [b_FF(1), b_FC(1)];    
% end


%keyboard;

return