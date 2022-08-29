home;
close all;
clear all;

%%
%define some colors
gold = [255,215,0]/255;
orange = [255, 165, 0]/255;
grey = [255, 255, 255]/2/255;

%% load data
load('exp4.mat')

[vmr_sign, idx, num_trials, tgt_dir, num_tgt, xp_all, yp_all, xv_all, yv_all, MV, Hang, Hdist,...
    Hstop, Hend, RT, Tang, Vang, CC, EC, Sess, md] = deal(exp4.vmr_sign, exp4.idx, exp4.num_trials, exp4.tgt_dir,...
    exp4.num_tgt, exp4.xp_all, exp4.yp_all, exp4.xv_all, exp4.yv_all, exp4.MV, exp4.Hang,...
    exp4.Hdist, exp4.Hstop, exp4.Hend, exp4.RT, exp4.Tang, exp4.Vang, exp4.CC, exp4.EC,...
    exp4.Sess, exp4.md);
nsubj = length(vmr_sign);


%% calculate adaptive response and learning rates

%NOTE: Null trials vs baseline trials give essentially the same result -- so use bln data
% figure; hold on; xlabel('Null movement directions (deg)'); ylabel('All baseline movement directions (deg)');
% plot(md_null_avg, md_bln_avg, 'b.', 'markersize', 10);

%calculate the mean baseline angle
md_null_avg = cellfun(@nanmean, md.null.all, 'UniformOutput', false);
md_bln_avg = cellfun(@nanmean, md.bln.all,'UniformOutput', false);

%optionally subtract baseline data (doesnt affect learning rate unless bias changes over trials)
bln_sub_flag = 0;
if bln_sub_flag ==1
    for i=1:nsubj
        for k=1:num_tgt
            md.CC.post.P(:,i,k) = md.CC.post.P(:,i,k)-md_bln_avg{k}(i);
            md.CC.post.N(:,i,k) = md.CC.post.N(:,i,k)-md_bln_avg{k}(i);
            
            md.CC.pre.P(:,i,k) = md.CC.pre.P(:,i,k)-md_bln_avg{k}(i);
            md.CC.pre.N(:,i,k) = md.CC.pre.N(:,i,k)-md_bln_avg{k}(i);
            
            md.CC.pert.P(:,i,k) = md.CC.pert.P(:,i,k)-md_bln_avg{k}(i);
            md.CC.pert.N(:,i,k) = md.CC.pert.N(:,i,k)-md_bln_avg{k}(i);
        end
        
    end
else
    md.CC.post.P = md.CC.post.P;
    md.CC.post.N = md.CC.post.N;
    
    md.CC.pre.N = md.CC.pre.N;
    md.CC.pre.P = md.CC.pre.P;
    
    md.CC.pert.N = md.CC.pert.N;
    md.CC.pert.P = md.CC.pert.P;
end

%calculate the AR
ar.CC.P = md.CC.post.P - md.CC.pert.P;
ar.CC.N = md.CC.post.N - md.CC.pert.N;
ar.CC.comb = ar.CC.P/2 + -ar.CC.N/2;
ar.CC.all = [ar.CC.P; ar.CC.N];

%calculate the AR across trials
ar_sub_avg.CC.P = squeeze(nanmean(ar.CC.P,1));
ar_sub_avg.CC.N = squeeze(nanmean(ar.CC.N,1));
ar_sub_avg.CC.comb = squeeze(nanmean(ar.CC.comb,1));
ar_sub_avg.CC.all = squeeze(nanmean(ar.CC.all,1));

%calculate learning rate on individual trials
LR.CC.P = ar.CC.P / 7.5;
LR.CC.N = ar.CC.N / -7.5;
LR.CC.comb = LR.CC.P/2 + LR.CC.N/2;
LR.CC.all = [LR.CC.P; LR.CC.N];

%calculate learning rate across trials
LR_sub_avg.CC.P = squeeze(nanmean(LR.CC.P,1));
LR_sub_avg.CC.N = squeeze(nanmean(LR.CC.N,1));
LR_sub_avg.CC.comb = squeeze(nanmean(LR.CC.comb,1));
LR_sub_avg.CC.all = squeeze(nanmean(LR.CC.all,1));

%calculate popluation mean and SE across subjects for each target direction
LR_gm.CC.P = mean(LR_sub_avg.CC.P,1);
LR_gm.CC.N = mean(LR_sub_avg.CC.N,1);
LR_gm.CC.comb = mean(LR_sub_avg.CC.comb,1);
LR_gm.CC.all = mean(LR_sub_avg.CC.all,1);

LR_se.CC.P = std(LR_sub_avg.CC.P,1)/sqrt(nsubj);
LR_se.CC.N = std(LR_sub_avg.CC.N,1)/sqrt(nsubj);
LR_se.CC.comb = std(LR_sub_avg.CC.comb,1)/sqrt(nsubj);
LR_se.CC.all = std(LR_sub_avg.CC.all,1)/sqrt(nsubj);

%% plot panel c: a bar plot showing the mean data and the 90 degree data

XL = {'90 degree data', 'Average across all target dir'};
YL = ['Learning rate'];
FT = ['Learning rate for different target directions (Takuji data)'];
color_order = repmat(ones(1,3)/2, num_tgt, 1);
d90_idx = find(tgt_dir==90);

paper_data_avg = [mean(mean(LR_sub_avg.CC.comb(:, d90_idx),2)), mean(mean(LR_sub_avg.CC.comb,2))];
paper_data_se = [std(mean(LR_sub_avg.CC.comb(:, d90_idx),2))/sqrt(nsubj), std(mean(LR_sub_avg.CC.comb,2))/sqrt(nsubj)];

h2 = display_LR_summary(paper_data_avg, paper_data_se, color_order, XL, YL, FT, []);

%get p value
[~, p_comp] = ttest(mean(LR_sub_avg.CC.comb(:, d90_idx),2), mean(LR_sub_avg.CC.comb,2));

%% plot panel b: a polar plot for all tgt directions
figure;
hold off;
%tgt_dir_W = wrapTo360([tgt_dir; tgt_dir(1)]);
tgt_dir_W = ([tgt_dir; tgt_dir(1)+360]);
ar_gm.CC.comb_W = [mean(ar_sub_avg.CC.comb,1), mean(ar_sub_avg.CC.comb(:,1))];
ar_se.CC.comb_W = [std(ar_sub_avg.CC.comb,0,1)/sqrt(nsubj), std(ar_sub_avg.CC.comb(:,1))/sqrt(nsubj)];

polarplot_interp(tgt_dir_W*pi/180,ar_gm.CC.comb_W,'m','linewidth',2); hold on;
polarplot_interp(tgt_dir_W*pi/180,ar_gm.CC.comb_W+ar_se.CC.comb_W,'m','linewidth',2, 'linestyle', '--'); hold on;
polarplot_interp(tgt_dir_W*pi/180,ar_gm.CC.comb_W-ar_se.CC.comb_W,'m','linewidth',2, 'linestyle', '--'); hold on;
polarplot(tgt_dir_W*pi/180,ar_gm.CC.comb_W,'.k','markersize',20); hold on;
%standard_error_shading_07_16_2015(ar_gm.CC.comb_W, ar_se.CC.comb_W, tgt_dir_W*pi/180, 1, 'm');
%polarplot_interp(tgt_dir_W*pi/180,tgt_dir_W*0+3.75,'k--','markersize',20); hold on;

%restrict to an adaptive response corresponding to 80% of full
rlim_plot = 0.8*7.5;

ax = gca;
ax.RTick = [0, 0.2, 0.4, 0.6, 0.8]*7.5;
ax.RTickLabel = '';
%ax.RGrid = 'off';
%ax.ThetaGrid = 'off';
ax.ThetaTick = sort(wrapTo360(tgt_dir));
rlim([0, rlim_plot])

title('Polar plot for all data');
shg


