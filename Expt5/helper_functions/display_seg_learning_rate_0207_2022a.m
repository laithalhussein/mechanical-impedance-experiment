function [ar_avg_sub, ar_sd_sub, ar_avg_gm, ar_avg_se, LR_sub, LR_gm, Rsq] =...
    display_seg_learning_rate_0207_2022a(AR, sub_id, title_prefix, bias_ar_flag)

V_err = [-7.5, -3.75, 0, 3.75, 7.5];
C_err = [-7.5, -3.75, 3.75, 7.5];

purple = [0.5,0,0.5];
grey = [0.5,0.5,0.5];
orange = [255, 165, 0]/255;
dgreen = [34, 139, 34]/255;
brown = [210,105,30]/255;
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
yellow = [255, 179, 0]/255;
dark_purple = [48, 25, 52]/255;
light_purple = [203, 195, 227]/255;
%% calculate the mean and SD of the adaptive response for each field

num_subjects = length(sub_id);

%start with NH
ar_avg.NH.V.P0.all = nanmean(AR.NH.V.P0.all,1);
ar_avg.NH.V.P3.P = nanmean(AR.NH.V.P3.P,1);
ar_avg.NH.V.P3.N = nanmean(AR.NH.V.P3.N,1);
ar_avg.NH.V.P7.P = nanmean(AR.NH.V.P7.P,1);
ar_avg.NH.V.P7.N = nanmean(AR.NH.V.P7.N,1);

ar_sd.NH.V.P0.all = nanstd(AR.NH.V.P0.all,0,1);
ar_sd.NH.V.P3.P = nanstd(AR.NH.V.P3.P,0,1);
ar_sd.NH.V.P3.N = nanstd(AR.NH.V.P3.N,0,1);
ar_sd.NH.V.P7.P = nanstd(AR.NH.V.P7.P,0,1);
ar_sd.NH.V.P7.N = nanstd(AR.NH.V.P7.N,0,1);

ar_avg.NH.C.P3.P = nanmean(AR.NH.C.P3.P,1);
ar_avg.NH.C.P3.N = nanmean(AR.NH.C.P3.N,1);
ar_avg.NH.C.P7.P = nanmean(AR.NH.C.P7.P,1);
ar_avg.NH.C.P7.N = nanmean(AR.NH.C.P7.N,1);

ar_sd.NH.C.P3.P = nanstd(AR.NH.C.P3.P,0,1);
ar_sd.NH.C.P3.N = nanstd(AR.NH.C.P3.N,0,1);
ar_sd.NH.C.P7.P = nanstd(AR.NH.C.P7.P,0,1);
ar_sd.NH.C.P7.N = nanstd(AR.NH.C.P7.N,0,1);

%do HI
ar_avg.HI.P3.P = nanmean(AR.HI.P3.P,1);
ar_avg.HI.P3.N = nanmean(AR.HI.P3.N,1);
ar_avg.HI.P7.P = nanmean(AR.HI.P7.P,1);
ar_avg.HI.P7.N = nanmean(AR.HI.P7.N,1);

ar_sd.HI.P3.P = nanstd(AR.HI.P3.P,0,1);
ar_sd.HI.P3.N = nanstd(AR.HI.P3.N,0,1);
ar_sd.HI.P7.P = nanstd(AR.HI.P7.P,0,1);
ar_sd.HI.P7.N = nanstd(AR.HI.P7.N,0,1);

%do HD
ar_avg.HD.P3.P = nanmean(AR.HD.P3.P,1);
ar_avg.HD.P3.N = nanmean(AR.HD.P3.N,1);
ar_avg.HD.P7.P = nanmean(AR.HD.P7.P,1);
ar_avg.HD.P7.N = nanmean(AR.HD.P7.N,1);

ar_sd.HD.P3.P = nanstd(AR.HD.P3.P,0,1);
ar_sd.HD.P3.N = nanstd(AR.HD.P3.N,0,1);
ar_sd.HD.P7.P = nanstd(AR.HD.P7.P,0,1);
ar_sd.HD.P7.N = nanstd(AR.HD.P7.N,0,1);

%% concatinate each condition into a matrix

%NH
ar_avg_sub.NH.V = [ar_avg.NH.V.P7.N; ar_avg.NH.V.P3.N; ar_avg.NH.V.P0.all; ar_avg.NH.V.P3.P; ar_avg.NH.V.P7.P];
ar_sd_sub.NH.V = [ar_sd.NH.V.P7.N; ar_sd.NH.V.P3.N; ar_sd.NH.V.P0.all; ar_sd.NH.V.P3.P; ar_sd.NH.V.P7.P];
ar_avg_sub.NH.C = [ar_avg.NH.C.P7.N; ar_avg.NH.C.P3.N; ar_avg.NH.C.P3.P; ar_avg.NH.C.P7.P];
ar_sd_sub.NH.C = [ar_sd.NH.C.P7.N; ar_sd.NH.C.P3.N; ar_sd.NH.C.P3.P; ar_sd.NH.C.P7.P];

%HI
% ar_avg_sub.HI = [ar_avg.HI.P7.N; ar_avg.HI.P3.N; ar_avg.NH.V.P0.all; ar_avg.HI.P3.P; ar_avg.HI.P7.P];
% ar_sd_sub.HI = [ar_sd.HI.P7.N; ar_sd.HI.P3.N; ar_sd.NH.V.P0.all; ar_sd.HI.P3.P; ar_sd.HI.P7.P];
ar_avg_sub.HI = [ar_avg.HI.P7.N; ar_avg.HI.P3.N; ar_avg.HI.P3.P; ar_avg.HI.P7.P];
ar_sd_sub.HI = [ar_sd.HI.P7.N; ar_sd.HI.P3.N; ar_sd.HI.P3.P; ar_sd.HI.P7.P];

%HD
% ar_avg_sub.HD = [ar_avg.HD.P7.N; ar_avg.HD.P3.N; ar_avg.NH.V.P0.all; ar_avg.HD.P3.P; ar_avg.HD.P7.P];
% ar_sd_sub.HD = [ar_sd.HD.P7.N; ar_sd.HD.P3.N; ar_sd.NH.V.P0.all; ar_sd.HD.P3.P; ar_sd.HD.P7.P];
ar_avg_sub.HD = [ar_avg.HD.P7.N; ar_avg.HD.P3.N; ar_avg.HD.P3.P; ar_avg.HD.P7.P];
ar_sd_sub.HD = [ar_sd.HD.P7.N; ar_sd.HD.P3.N; ar_sd.HD.P3.P; ar_sd.HD.P7.P];

%% optionally subtract the bias for each trial type independently
if bias_ar_flag
    ar_avg_sub.NH.V = ar_avg_sub.NH.V - mean(nanmean(ar_avg_sub.NH.V));
    ar_avg_sub.NH.C = ar_avg_sub.NH.C - mean(nanmean(ar_avg_sub.NH.C));
    ar_avg_sub.HI = ar_avg_sub.HI - mean(nanmean(ar_avg_sub.HI));
    ar_avg_sub.HD = ar_avg_sub.HD - mean(nanmean(ar_avg_sub.HD));    
end

%% calculate population avergae and standard error across subjects

ar_avg_gm.NH.V = mean(ar_avg_sub.NH.V,2);
ar_avg_se.NH.V = std(ar_avg_sub.NH.V,0,2)/sqrt(num_subjects);

ar_avg_gm.NH.C = mean(ar_avg_sub.NH.C,2);
ar_avg_se.NH.C = std(ar_avg_sub.NH.C,0,2)/sqrt(num_subjects);

ar_avg_gm.HI = mean(ar_avg_sub.HI,2);
ar_avg_se.HI = std(ar_avg_sub.HI,0,2)/sqrt(num_subjects);

ar_avg_gm.HD = mean(ar_avg_sub.HD,2);
ar_avg_se.HD = std(ar_avg_sub.HD,0,2)/sqrt(num_subjects);

%% calculate the slope for each subject and for the overall population

V_err_reg = [V_err(:), V_err(:)*0+1];
C_err_reg = [C_err(:), C_err(:)*0+1];

%begin with individual subjects
for k=1:num_subjects
    btmp = regress(ar_avg_sub.NH.V(:,k), V_err_reg);
    LR_sub.NH.V(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.NH.C(:,k), C_err_reg);
    LR_sub.NH.C(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.HI(:,k), C_err_reg);
    LR_sub.HI(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.HD(:,k), C_err_reg);
    LR_sub.HD(:,k) = btmp;
end

%get learning rate based on population averages
[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.NH.V, V_err_reg);
LR_gm.NH.V = btmp;
Rsq.NH.V = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.NH.C, C_err_reg);
LR_gm.NH.C = btmp;
Rsq.NH.C = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.HI, C_err_reg);
LR_gm.HI = btmp;
Rsq.HI = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.HD, C_err_reg);
LR_gm.HD = btmp;
Rsq.HD = stmp(1);
%% plot the results!

figure; hold on;

%first make the legend
msize = 20;
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', 'r', 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', 'b', 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', light_purple, 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', dark_purple, 'Markersize', msize);

leg_tmp = legend({'VMR trials', 'Clamp trials',...
    'VMR \rightarrow VC trials',...
    'VC \rightarrow VMR trials'});
legend('Autoupdate', 'Off');

%VMR adaptation
plot(V_err, ar_avg_gm.NH.V, 'Linestyle', 'None', 'Marker', '.', 'color', 'r', 'Markersize', msize);
display_xy_error(V_err, ar_avg_gm.NH.V, [], ar_avg_se.NH.V, 'color', 'r');

%Clamp adaptation
plot(C_err, ar_avg_gm.NH.C, 'Linestyle', 'None', 'Marker', '.', 'color', 'b', 'Markersize', msize);
display_xy_error(C_err, ar_avg_gm.NH.C, [], ar_avg_se.NH.C, 'color', 'b');

%VMR -> VC adaptation
plot(C_err, ar_avg_gm.HI, 'Linestyle', 'None', 'Marker', '.', 'color', light_purple, 'Markersize', msize);
display_xy_error(C_err, ar_avg_gm.HI, [], ar_avg_se.HI, 'color', light_purple);

%VC -> VMR adaptation
plot(C_err, ar_avg_gm.HD, 'Linestyle', 'None', 'Marker', '.', 'color', dark_purple, 'Markersize', msize);
display_xy_error(C_err, ar_avg_gm.HD, [], ar_avg_se.HD, 'color', dark_purple);

%Draw the regression lines for each case
plot(V_err, V_err_reg*LR_gm.NH.V, 'Linestyle', '-', 'color', 'r', 'Linewidth', 1);
plot(C_err, C_err_reg*LR_gm.NH.C, 'Linestyle', '-', 'color', 'b', 'Linewidth', 1);
plot(C_err, C_err_reg*LR_gm.HI, 'Linestyle', '-', 'color', light_purple, 'Linewidth', 1);
plot(C_err, C_err_reg*LR_gm.HD, 'Linestyle', '-', 'color', dark_purple, 'Linewidth', 1);

%plot y = -x
%plot(V_err, -V_err, 'Linestyle', '--', 'Color', 'k');

set(gca, 'xTick', V_err);
%grid on;
ax = gca;
ax.YTick = ax.XTick;
ylim([-6,6]);

title([title_prefix, ' learning']);
ylabel('Adaptive response (deg)');
xlabel('Errors (deg)');

set(leg_tmp, 'Location','Best');
return