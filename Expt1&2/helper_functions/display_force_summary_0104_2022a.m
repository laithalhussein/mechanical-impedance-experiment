function [LR_reg_sub_all] = display_force_summary_0104_2022a(f, tt)

%displays the average difference in positive and negative force profiles for the perturbation and adaptive response
%Also returns the learning rate based on regression

grey = [0.5,0.5,0.5];
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
orange = [255, 165, 0]/255;
dgreen = [34, 139, 34]/255;
brown = [210,105,30]/255;
fgain = 0.9;
%fgain =1;

[~, num_subjects, num_samples] = size(f.pert.FF.all);

%we will save the learning rate measurements, and the subject averaged
%ARs and ideal force data
LR_reg_sub_all = nan(2, num_subjects);

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
figure; hold on;
plot(tt, fgain*mean(FF_pert_diff_sub), 'color', orange); %ideal
plot(tt, mean(FF_AR_diff_sub), 'color', 'r');
standard_error_shading_07_16_2015(fgain*mean(FF_pert_diff_sub),std(FF_pert_diff_sub),tt,num_subjects,orange);
standard_error_shading_07_16_2015(mean(FF_AR_diff_sub),std(FF_AR_diff_sub),tt,num_subjects,'r');
ylim([-1,5]);
xlabel('Time wrt Target acquisition (sec)');
ylabel('Force (N)');
title('FF perturbation and adaptive response');
ax = gca;
ax.YTick = [0, 5];
ax.XTick = [-0.4, -0.2, 0];

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
for nsub=1:num_subjects
    reg_input_FF = [FF_pert_diff_sub(nsub,:); FF_pert_diff_sub(nsub,:)*0+1]'; %regress the ideal onto the adaptive response
    reg_input_FC = [FC_pert_diff_sub(nsub,:); FC_pert_diff_sub(nsub,:)*0+1]';
    
    b_FF = regress(FF_AR_diff_sub(nsub,:)', reg_input_FF);
    b_FC = regress(FC_AR_diff_sub(nsub,:)', reg_input_FC);
    
%     if b_FF(1) < 0.05, b_FF(1) = 0.1; end
%     if b_FF(1) > 0.25, b_FF(1) = b_FF(1) - 0.1; end
%     if b_FC(1) > 1, b_FC(1) = 0.9; end
%     if b_FC(1) < 0.26, b_FC(1) = b_FC(1)*2; end
    
    LR_reg_sub_all(:, nsub) = [b_FF(1), b_FC(1)];    
end

%get peak forces
peak_FF = max(fgain*FF_pert_diff_sub, [], 2);
peak_FC = max(FC_pert_diff_sub, [], 2);

peak_FF_avg = mean(peak_FF);
peak_FF_se = std(peak_FF)/sqrt(num_subjects);

peak_FC_avg = mean(peak_FC);
peak_FC_se = std(peak_FC)/sqrt(num_subjects);

%alternate way of calculating the peak forces
peak_FF2_avg = max(nanmean(fgain*FF_pert_diff_sub,1));
peak_FC2_avg = max(nanmean(FC_pert_diff_sub,1));

%get correlation coefficient
r_FF_sub = nan(num_subjects, 1);
r_FC_sub = nan(num_subjects, 1);

for nsub=1:num_subjects
    [r_FF_sub(nsub), ~] = corr(FF_pert_diff_sub(nsub,:)', FF_AR_diff_sub(nsub,:)');
    [r_FC_sub(nsub), ~] = corr(FC_pert_diff_sub(nsub,:)', FC_AR_diff_sub(nsub,:)');
end

%get SE and p value based on individual subjects, but mean based on pop. avg
r_FF_gm = corr(mean(FF_pert_diff_sub)', mean(FF_AR_diff_sub)');
r_FC_gm = corr(mean(FC_pert_diff_sub)', mean(FC_AR_diff_sub)');

f_FF_se = std(r_FF_sub)/sqrt(num_subjects);
f_FC_se = std(r_FC_sub)/sqrt(num_subjects);

[~, p_FF_corr, ~, s_FF_corr] = ttest(r_FF_sub);
[~, p_FC_corr, ~, s_FC_corr] = ttest(r_FC_sub);

%keyboard;

return