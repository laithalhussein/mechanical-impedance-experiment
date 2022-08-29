clear all;
close all;
home;

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

%Flags/constants
iqr_filter_flag = 1;
iqr_num = 3;
bias_flag = 1;
filter_md_flag = 1;

%% load data
load('exp3.mat')

[xp, yp, xv, yv, mo_dist, rt, mt, mov_dir, explicit_dir, tv_max, tvt_max,...
    num_subjects, num_trials, sub_id, idx, xp_, yp_, xv_, yv_, mo_dist_, rt_,...
    mt_, mov_dir_, explicit_dir_, tv_max_, tvt_max_] = deal(exp3.xp, exp3.yp, exp3.xv, exp3.yv,...
    exp3.mo_dist, exp3.rt, exp3.mt, exp3.mov_dir, exp3.explicit_dir, exp3.tv_max, exp3.tvt_max,...
    exp3.num_subjects, exp3.num_trials, exp3.sub_id, exp3.idx, exp3.xp_, exp3.yp_, exp3.xv_, exp3.yv_, exp3.mo_dist_, exp3.rt_,...
    exp3.mt_, exp3.mov_dir_, exp3.explicit_dir_, exp3.tv_max_, exp3.tvt_max_);

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

%% Produce figures that show implicit and explicit adaptive responses (panel e)

[imp.ar_avg_sub, imp.ar_sd_sub, imp.ar_avg_gm,...
    imp.ar_avg_se, imp.regc_sub, imp.regc_gm, imp.Rsq] = display_learning_rate_0206_2022a(IA, sub_id, 'Implicit', 0);

[exp.ar_avg_sub, exp.ar_sd_sub, exp.ar_avg_gm,...
    exp.ar_avg_se, exp.regc_sub, exp.regc_gm, exp.Rsq] = display_learning_rate_0206_2022a(EA, sub_id, 'Explicit', 0);

[total.ar_avg_sub, total.ar_sd_sub, total.ar_avg_gm,...
    total.ar_avg_se, total.regc_sub, total.regc_gm, total.Rsq] = display_learning_rate_0206_2022a(TA, sub_id, 'Total', 0);
close(3);


%% show mean and SE of learning rates in a bar plot

imp_LR_avg = [mean(imp.regc_sub.LIE(1,:)), mean(imp.regc_sub.HIE(1,:)), mean(imp.regc_sub.MIX.V(1,:)), mean(imp.regc_sub.MIX.C(1,:))];
imp_LR_se = [std(imp.regc_sub.LIE(1,:)), std(imp.regc_sub.HIE(1,:)), std(imp.regc_sub.MIX.V(1,:)), std(imp.regc_sub.MIX.C(1,:))]...
    /sqrt(num_subjects);

color_order = [red; blue; orange; cyan];

XL = {'VMR trials from a LIE','Clamp trials from a HIE', 'VMR trials from a Mixed Env', 'Clamp trials from a Mixed Env'};
YL = ['Learning rate (flipped)'];
FT = ['Effects of impedance on implicit learning rate'];
h = display_LR_summary(-imp_LR_avg, imp_LR_se, color_order, XL, YL, FT, []);

%optionally we plot the data on top of the bar plot (and connect them with lines)
data_bar_flag = 0;
if data_bar_flag
    ax = gca;
    
    %get the x-positions on the bar plots
    xpos_tmp = ax.XTick;
    
    %concatinate the individual learning rate data
    y_tmp = -[imp.regc_sub.LIE(1,:); imp.regc_sub.HIE(1,:); imp.regc_sub.MIX.V(1,:); imp.regc_sub.MIX.C(1,:)]';
    
    %loop through and plot (and connect) each coniditon on the bar plot
    for k = 2:4
        plot_data_on_bar(xpos_tmp(k-1:k), y_tmp(:,k-1:k), 'marker', 'o', 'Linestyle', '--', 'MarkerFaceColor', 'k', 'Color', grey);
    end
    
    ylim([0, 0.8]);
end

%% make a scatter plot of the individual subject learning rates
%We plot both the individual envs and mixed cases on the same plot
%Also calculate the % individual subject difference for both conditions

figure; hold on;
xlabel('Learning rate on VMR trials');
ylabel('Learning rate on VC trials');

plot(nan, nan, 'marker', 'o', 'color', 'k', 'displayname', 'Blocked data', 'Linestyle', 'None');
plot(nan, nan, 'marker', 'o', 'color', grey, 'displayname', 'Intermixed data', 'Linestyle', 'None');
leg_tmp = legend('show');
legend('Autoupdate', 'Off');
set(leg_tmp, 'Location', 'Best');

%plot each datapoint separately
plot(-imp.regc_sub.LIE(1,:), -imp.regc_sub.HIE(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'k');
plot(-imp.regc_sub.MIX.V(1,:), -imp.regc_sub.MIX.C(1,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', grey);

%connect each condition within participants
xcond = [-imp.regc_sub.LIE(1,:); -imp.regc_sub.MIX.V(1,:)];
ycond = [-imp.regc_sub.HIE(1,:); -imp.regc_sub.MIX.C(1,:)];
plot(xcond, ycond, 'Linestyle', '--', 'color', light_grey);

%plot the grand means
plot(mean(-imp.regc_sub.LIE(1,:)), mean(-imp.regc_sub.HIE(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(mean(-imp.regc_sub.MIX.V(1,:)), mean(-imp.regc_sub.MIX.C(1,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', grey, 'MarkerEdgeColor', grey);

%connect the means
xavg_cond = mean(xcond,2);
yavg_cond = mean(ycond,2);
plot(xavg_cond, yavg_cond, 'Linestyle', '--', 'color', grey);

%show error elipse
%plot_error_ellipse(-imp.regc_sub.LIE(1,:),-imp.regc_sub.HIE(1,:),black,1.5, 'show_CI');
%plot_error_ellipse(-imp.regc_sub.MIX.V(1,:),-imp.regc_sub.MIX.C(1,:),grey,1.5, 'show_CI');

%plot y = x
plot([0,1], [0,1], 'Linestyle', '--', 'Color', 'k');

%get correlation coefficient
[r.MIX, p.MIX] = corr(-imp.regc_sub.LIE(1,:)', -imp.regc_sub.HIE(1,:)');
[r.ISO, p.ISO] = corr(-imp.regc_sub.MIX.V(1,:)', -imp.regc_sub.MIX.C(1,:)', 'rows', 'complete');

set(leg_tmp, 'Location', 'Best');
xlim([0,0.8]);
ylim([0,0.8]);
title('Individual subject learning rates for all conditions');
%title(['Relationship between HIE and LIE learning rates across subjects, r = ', num2str(corr]);

ax = gca;
ax.XTick = [0, 0.2, 0.4, 0.6, 0.8];
ax.YTick = [0, 0.2, 0.4, 0.6, 0.8];

%% calculate fraction of implicit learning to total

%do this by regression
% reg_tmp = [imp.ar_avg_gm.MIX.V, imp.ar_avg_gm.MIX.V*0+1];
%  [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.MIX.V, reg_tmp);
% frac.MIX.V = stmp(1);

reg_tmp = [imp.ar_avg_gm.MIX.C, imp.ar_avg_gm.MIX.C*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.MIX.C, reg_tmp);
frac.MIX.C = stmp(1);

reg_tmp = [imp.ar_avg_gm.LIE, imp.ar_avg_gm.LIE*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.LIE, reg_tmp);
frac.LIE = stmp(1);

reg_tmp = [imp.ar_avg_gm.HIE, imp.ar_avg_gm.HIE*0+1];
 [~, ~, ~, ~, stmp]= regress(total.ar_avg_gm.HIE, reg_tmp);
frac.HIE = stmp(1);

%% print p-values
%%% p-values for seeing if learning rates are significantly less than 0 (its negative here)
[~, imp_pz.LIE, ~, imp_tstat.LIE] = ttest(imp.regc_sub.LIE(1,:), 0, 'tail', 'both');
[~, imp_pz.MIX.V, ~, imp_tstat.MIX.V] = ttest(imp.regc_sub.MIX.V(1,:), 0, 'tail', 'both');
[~, imp_pz.MIX.C, ~, imp_tstat.MIX.C] = ttest(imp.regc_sub.MIX.C(1,:), 0, 'tail', 'both');
[~, imp_pz.HIE, ~, imp_tstat.HIE] = ttest(imp.regc_sub.HIE(1,:), 0, 'tail', 'both');

%%% p-values for pairwise difference between HIE vs LIE and MIX.C vs MIX.V
[~, imp_pdiff1, ~, imp_tstat_diff1] = ttest(imp.regc_sub.LIE(1,:), imp.regc_sub.HIE(1,:));
[~, imp_pdiff2, ~, imp_tstat_diff2] = ttest(imp.regc_sub.MIX.V(1,:), imp.regc_sub.MIX.C(1,:));

%print the p-values
T = table(imp_pz.LIE, imp_pz.MIX.V, imp_pz.MIX.C, imp_pz.HIE, imp_pdiff1, imp_pdiff2);
T.Properties.VariableNames = {'imp_pz_LIE','imp_pz_MIXV', 'imp_pz_MIXC', 'imp_pz_HIE', 'HIE_vs_LIE', 'MIXV_vs_MIXC'};
disp(T)

%%%%repeat for explicit learning
exp_LR_avg = [mean(exp.regc_sub.LIE(1,:)), mean(exp.regc_sub.HIE(1,:)), mean(exp.regc_sub.MIX.V(1,:)), mean(exp.regc_sub.MIX.C(1,:))];
exp_LR_se = [std(exp.regc_sub.LIE(1,:)), std(exp.regc_sub.HIE(1,:)), std(exp.regc_sub.MIX.V(1,:)), std(exp.regc_sub.MIX.C(1,:))]...
    /sqrt(num_subjects);

FT2 = ['Effects of impedance on explicit learning rate'];
h = display_LR_summary(-exp_LR_avg, exp_LR_se, color_order, XL, YL, FT2, []);

%%% p-values for seeing if learning rates are significantly different than 0
[~, exp_pz.LIE, ~, exp_tstat.LIE] = ttest(exp.regc_sub.LIE(1,:), 0, 'tail', 'both');
[~, exp_pz.MIX.V, ~, exp_tstat.MIX.V] = ttest(exp.regc_sub.MIX.V(1,:), 0, 'tail', 'both');
[~, exp_pz.MIX.C, ~, exp_tstat.MIX.C] = ttest(exp.regc_sub.MIX.C(1,:), 0, 'tail', 'both');
[~, exp_pz.HIE, ~, exp_tstat.HIE] = ttest(exp.regc_sub.HIE(2,:), 0, 'tail', 'both');

%%% p-values for pairwise difference between HIE vs LIE and MIX.C vs MIX.V
[~, exp_pdiff1, ~, exp_tstat_diff1] = ttest(exp.regc_sub.LIE(1,:), exp.regc_sub.HIE(1,:));
[~, exp_pdiff2, ~, exp_tstat_diff2] = ttest(exp.regc_sub.MIX.V(1,:), exp.regc_sub.MIX.C(1,:));

%print the p-values in a table
