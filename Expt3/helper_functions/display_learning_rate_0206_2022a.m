function [ar_avg_sub, ar_sd_sub, ar_avg_gm, ar_avg_se, LR_sub, LR_gm, Rsq] =...
    display_learning_rate_0206_2022a(AR, sub_id, title_prefix, bias_ar_flag)

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
%% calculate the mean and SD of the adaptive response for each field

num_subjects = length(sub_id);

%start with LIE
ar_avg.LIE.V.P3.P = nanmean(AR.LIE.V.P3.P,1);
ar_avg.LIE.V.P3.N = nanmean(AR.LIE.V.P3.N,1);
ar_avg.LIE.V.P7.P = nanmean(AR.LIE.V.P7.P,1);
ar_avg.LIE.V.P7.N = nanmean(AR.LIE.V.P7.N,1);
ar_avg.LIE.V.P0.all = nanmean(AR.LIE.V.P0.all,1);

ar_sd.LIE.V.P3.P = nanstd(AR.LIE.V.P3.P,0,1);
ar_sd.LIE.V.P3.N = nanstd(AR.LIE.V.P3.N,0,1);
ar_sd.LIE.V.P7.P = nanstd(AR.LIE.V.P7.P,0,1);
ar_sd.LIE.V.P7.N = nanstd(AR.LIE.V.P7.N,0,1);
ar_sd.LIE.V.P0.all = nanstd(AR.LIE.V.P0.all,0,1);

%do MIXED
ar_avg.MIX.V.P3.P = nanmean(AR.MIX.V.P3.P,1);
ar_avg.MIX.V.P3.N = nanmean(AR.MIX.V.P3.N,1);
ar_avg.MIX.V.P7.P = nanmean(AR.MIX.V.P7.P,1);
ar_avg.MIX.V.P7.N = nanmean(AR.MIX.V.P7.N,1);
ar_avg.MIX.V.P0.all = nanmean(AR.MIX.V.P0.all,1);

ar_sd.MIX.V.P3.P = nanstd(AR.MIX.V.P3.P,0,1);
ar_sd.MIX.V.P3.N = nanstd(AR.MIX.V.P3.N,0,1);
ar_sd.MIX.V.P7.P = nanstd(AR.MIX.V.P7.P,0,1);
ar_sd.MIX.V.P7.N = nanstd(AR.MIX.V.P7.N,0,1);
ar_sd.MIX.V.P0.all = nanstd(AR.MIX.V.P0.all,0,1);

ar_avg.MIX.C.P3.P = nanmean(AR.MIX.C.P3.P,1);
ar_avg.MIX.C.P3.N = nanmean(AR.MIX.C.P3.N,1);
ar_avg.MIX.C.P7.P = nanmean(AR.MIX.C.P7.P,1);
ar_avg.MIX.C.P7.N = nanmean(AR.MIX.C.P7.N,1);

ar_sd.MIX.C.P3.P = nanstd(AR.MIX.C.P3.P,0,1);
ar_sd.MIX.C.P3.N = nanstd(AR.MIX.C.P3.N,0,1);
ar_sd.MIX.C.P7.P = nanstd(AR.MIX.C.P7.P,0,1);
ar_sd.MIX.C.P7.N = nanstd(AR.MIX.C.P7.N,0,1);
    
%do HIE
ar_avg.HIE.C.P3.P = nanmean(AR.HIE.C.P3.P,1);
ar_avg.HIE.C.P3.N = nanmean(AR.HIE.C.P3.N,1);
ar_avg.HIE.C.P7.P = nanmean(AR.HIE.C.P7.P,1);
ar_avg.HIE.C.P7.N = nanmean(AR.HIE.C.P7.N,1);
ar_avg.HIE.V.P0.all = nanmean(AR.HIE.V.P0.all,1);

ar_sd.HIE.C.P3.P = nanstd(AR.HIE.C.P3.P,0,1);
ar_sd.HIE.C.P3.N = nanstd(AR.HIE.C.P3.N,0,1);
ar_sd.HIE.C.P7.P = nanstd(AR.HIE.C.P7.P,0,1);
ar_sd.HIE.C.P7.N = nanstd(AR.HIE.C.P7.N,0,1);
ar_sd.HIE.V.P0.all = nanstd(AR.HIE.V.P0.all,0,1);

%% concatinate each condition into a matrix

% ar_avg_sub.LIE = [ar_avg.LIE.V.P7.N; ar_avg.LIE.V.P3.N; ar_avg.LIE.V.P0.all; ar_avg.LIE.V.P3.P; ar_avg.LIE.V.P7.P];
% ar_sd_sub.LIE = [ar_sd.LIE.V.P7.N; ar_sd.LIE.V.P3.N; ar_sd.LIE.V.P0.all; ar_sd.LIE.V.P3.P; ar_sd.LIE.V.P7.P];
ar_avg_sub.LIE = [ar_avg.LIE.V.P7.N; ar_avg.LIE.V.P3.N; ar_avg.LIE.V.P3.P; ar_avg.LIE.V.P7.P];
ar_sd_sub.LIE = [ar_sd.LIE.V.P7.N; ar_sd.LIE.V.P3.N; ar_sd.LIE.V.P3.P; ar_sd.LIE.V.P7.P];

% ar_avg_sub.MIX.V = [ar_avg.MIX.V.P7.N; ar_avg.MIX.V.P3.N; ar_avg.MIX.V.P0.all; ar_avg.MIX.V.P3.P; ar_avg.MIX.V.P7.P];
% ar_sd_sub.MIX.V = [ar_sd.MIX.V.P7.N; ar_sd.MIX.V.P3.N; ar_sd.MIX.V.P0.all; ar_sd.MIX.V.P3.P; ar_sd.MIX.V.P7.P];
ar_avg_sub.MIX.V = [ar_avg.MIX.V.P7.N; ar_avg.MIX.V.P3.N; ar_avg.MIX.V.P3.P; ar_avg.MIX.V.P7.P];
ar_sd_sub.MIX.V = [ar_sd.MIX.V.P7.N; ar_sd.MIX.V.P3.N; ar_sd.MIX.V.P3.P; ar_sd.MIX.V.P7.P];
ar_avg_sub.MIX.C = [ar_avg.MIX.C.P7.N; ar_avg.MIX.C.P3.N; ar_avg.MIX.C.P3.P; ar_avg.MIX.C.P7.P];
ar_sd_sub.MIX.C = [ar_sd.MIX.C.P7.N; ar_sd.MIX.C.P3.N; ar_sd.MIX.C.P3.P; ar_sd.MIX.C.P7.P];

% ar_avg_sub.HIE = [ar_avg.HIE.C.P7.N; ar_avg.HIE.C.P3.N; ar_avg.HIE.V.P0.all; ar_avg.HIE.C.P3.P; ar_avg.HIE.C.P7.P];
% ar_sd_sub.HIE = [ar_sd.HIE.C.P7.N; ar_sd.HIE.C.P3.N; ar_sd.HIE.V.P0.all; ar_sd.HIE.C.P3.P; ar_sd.HIE.C.P7.P];
ar_avg_sub.HIE = [ar_avg.HIE.C.P7.N; ar_avg.HIE.C.P3.N; ar_avg.HIE.C.P3.P; ar_avg.HIE.C.P7.P];
ar_sd_sub.HIE = [ar_sd.HIE.C.P7.N; ar_sd.HIE.C.P3.N; ar_sd.HIE.C.P3.P; ar_sd.HIE.C.P7.P];

%% optionally subtract the bias for each trial type independently
if bias_ar_flag
    ar_avg_sub.LIE = ar_avg_sub.LIE - mean(mean(ar_avg_sub.LIE));
    ar_avg_sub.HIE = ar_avg_sub.HIE - mean(mean(ar_avg_sub.HIE));
    ar_avg_sub.MIX.V = ar_avg_sub.MIX.V - mean(mean(ar_avg_sub.MIX.V));
    ar_avg_sub.MIX.C = ar_avg_sub.MIX.C - mean(mean(ar_avg_sub.MIX.C));    
end

%% calculate population avergae and standard error across subjects

ar_avg_gm.LIE = mean(ar_avg_sub.LIE,2);
ar_avg_se.LIE = std(ar_avg_sub.LIE,0,2)/sqrt(num_subjects);

ar_avg_gm.MIX.V = mean(ar_avg_sub.MIX.V,2);
ar_avg_se.MIX.V = std(ar_avg_sub.MIX.V,0,2)/sqrt(num_subjects);

ar_avg_gm.MIX.C = mean(ar_avg_sub.MIX.C,2);
ar_avg_se.MIX.C = std(ar_avg_sub.MIX.C,0,2)/sqrt(num_subjects);

ar_avg_gm.HIE = mean(ar_avg_sub.HIE,2);
ar_avg_se.HIE = std(ar_avg_sub.HIE,0,2)/sqrt(num_subjects);

%% calculate the slope for each subject and for the overall population

C_err_reg = [C_err(:), C_err(:)*0+1];
if size(ar_avg_sub.MIX.V,1)==length(C_err)
    V_err = C_err;    
end
V_err_reg = [V_err(:), V_err(:)*0+1];

%begin with individual subjects
for k=1:num_subjects
    btmp = regress(ar_avg_sub.LIE(:,k), V_err_reg);
    LR_sub.LIE(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.MIX.V(:,k), V_err_reg);
    LR_sub.MIX.V(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.MIX.C(:,k), C_err_reg);
    LR_sub.MIX.C(:,k) = btmp;
    
    btmp = regress(ar_avg_sub.HIE(:,k), C_err_reg);
    LR_sub.HIE(:,k) = btmp;
end

%get learning rate based on population averages
[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.LIE, V_err_reg);
LR_gm.LIE = btmp;
Rsq.LIE = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.MIX.V, V_err_reg);
LR_gm.MIX.V = btmp;
Rsq.MIX.V = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.MIX.C, C_err_reg);
LR_gm.MIX.C = btmp;
Rsq.MIX.C = stmp(1);

[btmp, ~, ~, ~, stmp] = regress(ar_avg_gm.HIE, C_err_reg);
LR_gm.HIE = btmp;
Rsq.HIE = stmp(1);

%% plot the results!

figure; hold on;

%first make the legend
msize = 20;
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', 'r', 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', 'b', 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', orange, 'Markersize', msize);
plot(NaN, NaN, 'Linestyle', 'None', 'Marker', '.', 'color', cyan, 'Markersize', msize);

leg_tmp = legend({'VMR perturbation from a LIE', 'Clamp perturbation from a HIE',...
    'VMR perturbations from a Mixed Env',...
    'Clamp perturbation from a Mixed Env'});
legend('Autoupdate', 'Off');

%LIE, VMR adaptation
plot(V_err, ar_avg_gm.LIE, 'Linestyle', 'None', 'Marker', '.', 'color', 'r', 'Markersize', msize);
display_xy_error(V_err, ar_avg_gm.LIE, [], ar_avg_se.LIE, 'color', 'r');

%HIE, clamp adaptation
plot(C_err, ar_avg_gm.HIE, 'Linestyle', 'None', 'Marker', '.', 'color', 'b', 'Markersize', msize);
display_xy_error(C_err, ar_avg_gm.HIE, [], ar_avg_se.HIE, 'color', 'b');

%MIX, VMR adaptation
plot(V_err, ar_avg_gm.MIX.V, 'Linestyle', 'None', 'Marker', '.', 'color', orange, 'Markersize', msize);
display_xy_error(V_err, ar_avg_gm.MIX.V, [], ar_avg_se.MIX.V, 'color', orange);

%MIX, clamp adaptation
plot(C_err, ar_avg_gm.MIX.C, 'Linestyle', 'None', 'Marker', '.', 'color', cyan, 'Markersize', msize);
display_xy_error(C_err, ar_avg_gm.MIX.C, [], ar_avg_se.MIX.C, 'color', cyan);

%Draw the regression lines for each case
plot(V_err, V_err_reg*LR_gm.LIE, 'Linestyle', '-', 'color', 'r', 'Linewidth', 1);
plot(C_err, C_err_reg*LR_gm.HIE, 'Linestyle', '-', 'color', 'b', 'Linewidth', 1);
plot(V_err, V_err_reg*LR_gm.MIX.V, 'Linestyle', '-', 'color', orange, 'Linewidth', 1);
plot(C_err, C_err_reg*LR_gm.MIX.C, 'Linestyle', '-', 'color', cyan, 'Linewidth', 1);

%plot y = -x
%plot(V_err, -V_err, 'Linestyle', '--', 'Color', 'k');

set(gca, 'xTick', V_err);
%grid on;
ax = gca;
ax.YTick = ax.XTick;
ylim([-5,5]);

title([title_prefix, ' learning']);
ylabel('Adaptive response (deg)');
xlabel('Errors (deg)');

set(leg_tmp, 'Location','Best');

return