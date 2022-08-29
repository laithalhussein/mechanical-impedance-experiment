clear all;
close all;
home;

load('HIE_seg_exp_dat_all_12_subjects.mat');
info = info_all;
dat = dat_all;

num_subjects = info.num_subjects;
num_trials = info.num_trials;

%Colors
blue = [0,0,1];
black = [0,0,0];
purple = [0.5,0,0.5];
grey = [0.5,0.5,0.5];
orange = [255, 165, 0]/255;
dgreen = [34, 139, 34]/255;
brown = [210,105,30]/255;
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
yellow = [255, 179, 0]/255;
red = [255, 0, 0]/255;
gold = [255,215,0]/255;
dark_purple = [48, 25, 52]/255;
light_purple = [203, 195, 227]/255;

%Flags/constants
iqr_filter_flag = 1;
iqr_num = 3;

%% load data

load('exp5.mat')

[xp, yp, xv, yv, mo_dist, rt, mt, mov_dir, explicit_dir, tv_max, tvt_max,...
    num_subjects, num_trials, sub_id, idx, xp_, yp_, xv_, yv_, mo_dist_, rt_,...
    mt_, mov_dir_, explicit_dir_, tv_max_, tvt_max_] = deal(exp5.xp, exp5.yp, exp5.xv, exp5.yv,...
    exp5.mo_dist, exp5.rt, exp5.mt, exp5.mov_dir, exp5.explicit_dir, exp5.tv_max, exp5.tvt_max,...
    exp5.num_subjects, exp5.num_trials, exp5.sub_id, exp5.idx, exp5.xp_, exp5.yp_, exp5.xv_, exp5.yv_, exp5.mo_dist_, exp5.rt_,...
    exp5.mt_, exp5.mov_dir_, exp5.explicit_dir_, exp5.tv_max_, exp5.tvt_max_);

%% Calculate adaptation

%calculate implicit movement direction
implicit_dir = mov_dir - explicit_dir;

%Calculate adaptation
Ndiff = 1; %Use 2 for post-pre, and 1 for post-pert
implicit_adaptation = [NaN(Ndiff,num_subjects); implicit_dir([1+Ndiff:end],:) - implicit_dir([1:end-Ndiff],:)];
explicit_adaptation = [NaN(Ndiff,num_subjects); explicit_dir([1+Ndiff:end],:) - explicit_dir([1:end-Ndiff],:)];
total_adaptation = [NaN(Ndiff,num_subjects); mov_dir([1+Ndiff:end],:) - mov_dir([1:end-Ndiff],:)];

%remove first and optionally second trial of each block
if Ndiff==1
    implicit_adaptation(idx.block_onset.all,:) = NaN;
    explicit_adaptation(idx.block_onset.all,:) = NaN;
    total_adaptation(idx.block_onset.all,:) = NaN;
elseif Ndiff==2
    implicit_adaptation([idx.block_onset.all, idx.block_onset.all+1],:) = NaN;
    explicit_adaptation([idx.block_onset.all, idx.block_onset.all+1],:) = NaN;
    total_adaptation([idx.block_onset.all, idx.block_onset.all+1],:) = NaN;
end

%To deal with the fact that we never observe an adaptive response to the error experienced on the very last trial,...
%we pad the data with an extra NaN
implicit_adaptation(end+1,:) = NaN;
explicit_adaptation(end+1,:) = NaN;
total_adaptation(end+1,:) = NaN;

%Segment the data based on the error that was experienced on the PREVIOUS trial
IA = extract_adaptation_idx(implicit_adaptation, idx);
EA = extract_adaptation_idx(explicit_adaptation, idx);
TA = extract_adaptation_idx(total_adaptation, idx);

% filter adaptation data

if iqr_filter_flag
    %Filter based on 
    IA = filter_adaptation_data(IA, num_subjects, iqr_num);
    EA = filter_adaptation_data(EA, num_subjects, iqr_num);
    TA = filter_adaptation_data(TA, num_subjects, iqr_num);    
end

%% Produce figures that show implicit and explicit adaptive responses (panel b)

[imp.ar_avg_sub, imp.ar_sd_sub, imp.ar_avg_gm,...
    imp.ar_avg_se, imp.regc_sub, imp.regc_gm, imp.Rsq] = display_seg_learning_rate_0207_2022a(IA, sub_id, 'Implicit', 1);

[exp.ar_avg_sub, exp.ar_sd_sub, exp.ar_avg_gm,...
    exp.ar_avg_se, exp.regc_sub, exp.regc_gm, exp.Rsq] = display_seg_learning_rate_0207_2022a(EA, sub_id, 'Explicit', 1);

[total.ar_avg_sub, total.ar_sd_sub, total.ar_avg_gm,...
    total.ar_avg_se, total.regc_sub, total.regc_gm, total.Rsq] = display_seg_learning_rate_0207_2022a(TA, sub_id, 'Total', 1);
close(3);

%% show mean and SE of learning rates in a bar plot (panel d)

imp_LR_avg = [mean(imp.regc_sub.NH.V(1,:)), mean(imp.regc_sub.NH.C(1,:)),...
    mean(imp.regc_sub.HI(1,:)), mean(imp.regc_sub.HD(1,:))];
imp_LR_se = [std(imp.regc_sub.NH.V(1,:)), std(imp.regc_sub.NH.C(1,:)), std(imp.regc_sub.HI(1,:)), std(imp.regc_sub.HD(1,:))]...
    /sqrt(num_subjects);

color_order = [red; blue; light_purple; dark_purple];

XL = {'VMR trials', 'Clamp trials', 'VMR \rightarrow VC trials', 'VC \rightarrow VMR trials'};
YL = ['Learning rate (flipped)'];
FT = ['Implicit learning rate for all trial types'];
h = display_LR_summary(-imp_LR_avg, imp_LR_se, color_order, XL, YL, FT, [-0.05, 0.6]);

%%% p-values for seeing if learning rates are significantly different than 0 (its negative here)
[~, imp_pz.NH.V, ~, imp_tstat.NH.V] = ttest(imp.regc_sub.NH.V(1,:), 0, 'tail', 'both');
[~, imp_pz.NH.C, ~, imp_tstat.NH.C] = ttest(imp.regc_sub.NH.C(1,:), 0, 'tail', 'both');
[~, imp_pz.HI, ~, imp_tstat.HI] = ttest(imp.regc_sub.HI(1,:), 0, 'tail', 'both');
[~, imp_pz.HD, ~, imp_tstat.HD] = ttest(imp.regc_sub.HD(1,:), 0, 'tail', 'both');

%%% p-values for pairwise difference between various trial types
[~, imp_pdiff1, ~, imp_tstat1] = ttest(imp.regc_sub.NH.V(1,:), imp.regc_sub.NH.C(1,:)); %VC vs VMR
[~, imp_pdiff2, ~, imp_tstat2] = ttest(imp.regc_sub.HI(1,:), imp.regc_sub.HD(1,:)); %VMR -> VC vs VC -> VMR
[~, imp_pdiff3, ~, imp_tstat3] = ttest(imp.regc_sub.NH.V(1,:), imp.regc_sub.HD(1,:)); %VMR vs VC -> VMR
[~, imp_pdiff4, ~, imp_tstat4] = ttest(imp.regc_sub.NH.C(1,:), imp.regc_sub.HI(1,:)); %VC vs VMR -> VC
[~, imp_pdiff5, ~, imp_tstat5] = ttest(imp.regc_sub.HI(1,:), imp.regc_sub.NH.V(1,:)); %VMR -> VC vs VMR
[~, imp_pdiff6, ~, imp_tstat6] = ttest(imp.regc_sub.HD(1,:), imp.regc_sub.NH.C(1,:)); %VC vs VC -> VMR

%print the p-values
T = table(imp_pz.NH.V, imp_pz.NH.C, imp_pz.HI, imp_pz.HD, imp_pdiff1, imp_pdiff2, imp_pdiff3, imp_pdiff4);
T.Properties.VariableNames = {'imp_pz_NHV', 'imp_pz_NHC', 'imp_pz_HI', 'imp_pz_HD', 'NHV_vs_NHC', 'HI_vs_HD',...
    'NHV_vs_HD', 'NHC_vs_HI'};
disp(T)

%%%%repeat for explicit learning
exp_LR_avg = [mean(exp.regc_sub.NH.V(1,:)), mean(exp.regc_sub.NH.C(1,:)),...
    mean(exp.regc_sub.HI(1,:)), mean(exp.regc_sub.HD(1,:))];
exp_LR_se = [std(exp.regc_sub.NH.V(1,:)), std(exp.regc_sub.NH.C(1,:)), std(exp.regc_sub.HI(1,:)), std(exp.regc_sub.HD(1,:))]...
    /sqrt(num_subjects);

FT2 = ['Explicit learning rate for all trial types'];
h = display_LR_summary(-exp_LR_avg, exp_LR_se, color_order, XL, YL, FT2, [-0.05, 0.6]);

%%% p-values for seeing if learning rates are significantly different from 0
[~, exp_pz.NH.V, ~, exp_tstat.NH.V] = ttest(exp.regc_sub.NH.V(1,:), 0, 'tail', 'both');
[~, exp_pz.NH.C, ~, exp_tstat.NH.C] = ttest(exp.regc_sub.NH.C(1,:), 0, 'tail', 'both');
[~, exp_pz.HI, ~, exp_tstat.HI] = ttest(exp.regc_sub.HI(1,:), 0, 'tail', 'both');
[~, exp_pz.HD, ~, exp_tstat.HD] = ttest(exp.regc_sub.HD(1,:), 0, 'tail', 'both');

%%% p-values for pairwise difference between various trial types
[~, exp_pdiff1, ~, exp_tstat1] = ttest(exp.regc_sub.NH.V(1,:), exp.regc_sub.NH.C(1,:)); %VC vs VMR
[~, exp_pdiff2, ~, exp_tstat2] = ttest(exp.regc_sub.HI(1,:), exp.regc_sub.HD(1,:)); %VMR -> VC vs VC -> VMR
[~, exp_pdiff3, ~, exp_tstat3] = ttest(exp.regc_sub.NH.V(1,:), exp.regc_sub.HD(1,:)); %VMR vs VC -> VMR
[~, exp_pdiff4, ~, exp_tstat4] = ttest(exp.regc_sub.NH.C(1,:), exp.regc_sub.HI(1,:)); %VC vs VMR -> VC
[~, exp_pdiff5, ~, exp_tstat5] = ttest(exp.regc_sub.HI(1,:), exp.regc_sub.NH.V(1,:)); %VMR -> VC vs VMR
[~, exp_pdiff6, ~, exp_tstat6] = ttest(exp.regc_sub.HD(1,:), exp.regc_sub.NH.C(1,:)); %VC vs VC -> VMR

%print the p-values
T = table(exp_pz.NH.V, exp_pz.NH.C, exp_pz.HI, exp_pz.HD);
T.Properties.VariableNames = {'exp_pz_NHV','exp_pz_NHC', 'exp_pz_HI', 'exp_pz_HD'};
disp(T)

%% calculate fraction of implicit learning to total

%do this by regression
reg_tmp = [imp.ar_avg_gm.NH.V, imp.ar_avg_gm.NH.V*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.NH.V, reg_tmp);
frac.NH.V = stmp(1);

reg_tmp = [imp.ar_avg_gm.NH.C, imp.ar_avg_gm.NH.C*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.NH.C, reg_tmp);
frac.NH.C = stmp(1);

reg_tmp = [imp.ar_avg_gm.HI, imp.ar_avg_gm.HI*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.HI, reg_tmp);
frac.HI = stmp(1);

reg_tmp = [imp.ar_avg_gm.HD, imp.ar_avg_gm.HD*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.HD, reg_tmp);
frac.HD = stmp(1);

%% determine and make a bar plot of contribution of first half vs second half of movement (panel e)

[b_seg_sub] = determine_segment_contribution2(-imp.regc_sub.NH.V(1,:), -imp.regc_sub.NH.C(1,:),...
    -imp.regc_sub.HD(1,:), -imp.regc_sub.HI(1,:));

b_seg_avg = mean(b_seg_sub,1);
b_seg_se = std(b_seg_sub,0,1)/sqrt(num_subjects);

%get p values
[~,p_seg,~,tstat_seg] = ttest(b_seg_sub(:,2), b_seg_sub(:,3));
[~,p_fh,~,tstat_fh] = ttest(b_seg_sub(:,2));

% show them in a bar plot
h = figure; hold on;
xb_spacing = 0.55;

num_conditions = 2;

xb1 = linspace(1,xb_spacing*num_conditions+xb_spacing,num_conditions);

for k=1:num_conditions
    barwitherr(NaN, xb1(k), b_seg_avg(k+1), xb_spacing-0.05, 'Edgecolor', 'k',...
        'FaceColor', color_order(k,:),'Linewidth', 2); hold on;
    display_xy_error(xb1(k), b_seg_avg(k+1), [], b_seg_se(k+1), 'color', 'k');
end

XL = {'Increase in learning due to presence of VC in first half', 'Increase in learning due to presence of VC in second half'};
xticks([xb1(:)]);
xticklabels(XL);
xtickangle(45)
set(gca,'layer','top','box','off');
ylabel('Increase in learning');
ylim([0, 0.3]);
%axis tight;

%grid on;
ax = gca;
% ax.YTick = [0, 0.5, 1];
ax.YTick = [0:0.1:0.5];

%% scatter plot of individual subject learning rates for trials that end with VMRs vs those that end in VCs (panel c)

figure; hold on;
xlabel('learning rate on full VMR trials or VC \rightarrow VMR trials');
ylabel('learning rate on full VC trials or VMR \rightarrow VC trials');

plot(nan, nan, 'marker', 'o', 'color', 'k', 'displayname', 'Full perturbation type', 'linestyle', 'none');
plot(nan, nan, 'marker', 'o', 'color', grey, 'displayname', 'Segmented perturbation type', 'linestyle', 'none');
leg_tmp = legend('show');
legend('Autoupdate', 'Off');
set(leg_tmp, 'Location', 'Best');

%plot each datapoint separately
plot(-imp.regc_sub.NH.V(1,:), -imp.regc_sub.NH.C(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'k');
plot(-imp.regc_sub.HD(1,:), -imp.regc_sub.HI(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', grey);

%connect each condition within participants
xcond = [-imp.regc_sub.NH.V(1,:); -imp.regc_sub.HD(1,:)]; %contanitating these as a matrix
ycond = [-imp.regc_sub.NH.C(1,:); -imp.regc_sub.HI(1,:)];
plot(xcond, ycond, 'Linestyle', '--', 'color', light_grey);

%plot the grand means
plot(mean(-imp.regc_sub.NH.V(1,:)), mean(-imp.regc_sub.NH.C(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(mean(-imp.regc_sub.HD(1,:)), mean(-imp.regc_sub.HI(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);

%connect te grand means
xavg_cond = mean(xcond,2);
yavg_cond = mean(ycond,2);
%plot(xavg_cond, yavg_cond, 'Linestyle', '--', 'color', light_grey);

%plot y = x
plot([0,1], [0,1], 'Linestyle', '--', 'Color', 'k');

%plot 95% CI Ellipses
% plot_error_ellipse(-imp.regc_sub.NH.V(1,:),-imp.regc_sub.NH.C(1,:),'k',1.5, 'show_CI');
% plot_error_ellipse(-imp.regc_sub.HD(1,:),-imp.regc_sub.HI(1,:),grey,1.5, 'show_CI');

xlim([0,0.8]);
ylim([0,0.8]);
title('Individual subject learning rates');
ax = gca;
ax.YTick = [0:0.2:0.8];
ax.XTick = [0:0.2:0.8];

%% also plot learning rates from perturbations that START with VMR against perturbations that START with VCs (panel c)

figure; hold on;
xlabel('Learning rate from perturbations starting in VMRs');
ylabel('Learning rate from perturbations starting in VCs');

plot(nan, nan, 'marker', 'o', 'color', 'k', 'displayname', 'Full perturbation type', 'linestyle', 'none');
plot(nan, nan, 'marker', 'o', 'color', grey, 'displayname', 'Segmented perturbation type', 'linestyle', 'none');
leg_tmp = legend('show');
legend('Autoupdate', 'Off');
set(leg_tmp, 'Location', 'Best');

%plot each datapoint separately
plot(-imp.regc_sub.NH.V(1,:), -imp.regc_sub.NH.C(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'k');
plot(-imp.regc_sub.HI(1,:), -imp.regc_sub.HD(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', grey);

%connect each condition within participants
xcond = [-imp.regc_sub.NH.V(1,:); -imp.regc_sub.HI(1,:)]; %contanitating these as a matrix
ycond = [-imp.regc_sub.NH.C(1,:); -imp.regc_sub.HD(1,:)];
plot(xcond, ycond, 'Linestyle', '--', 'color', light_grey);

%plot the grand means
plot(mean(-imp.regc_sub.NH.V(1,:)), mean(-imp.regc_sub.NH.C(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(mean(-imp.regc_sub.HI(1,:)), mean(-imp.regc_sub.HD(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);

%connect the grand means
xavg_cond = mean(xcond,2);
yavg_cond = mean(ycond,2);
%plot(xavg_cond, yavg_cond, 'Linestyle', '--', 'color', light_grey);

%plot y = x
plot([0,1], [0,1], 'Linestyle', '--', 'Color', 'k');

%plot 95% CI Ellipses
% plot_error_ellipse(-imp.regc_sub.NH.V(1,:),-imp.regc_sub.NH.C(1,:),'k',1.5, 'show_CI');
% plot_error_ellipse(-imp.regc_sub.HD(1,:),-imp.regc_sub.HI(1,:),grey,1.5, 'show_CI');

xlim([0,0.8]);
ylim([0,0.8]);
title('Individual subject learning rates');
ax = gca;
ax.YTick = [0:0.2:0.8];
ax.XTick = [0:0.2:0.8];