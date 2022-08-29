clear all;
close all;
home;

if ispc
    %cd('C:\Users\laith\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\mat_files\ALL_DATA'); %laptop
    cd('C:\Users\ryanm\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\mat_files\ALL_DATA'); %Ryan's computer!
else
    cd('/Users/laithalhussein/Dropbox (HNL)/Laith_files/HSE_exp/expb/3deg_data/mat_files/ALL_DATA'); %Macbook
end

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
if ispc
    hse_data=load(strcat(Cdr,'\',hse_filename)); %windows
else
    hse_data=load(strcat(Cdr,'/',hse_filename)); %mac
end
draw=hse_data.dat;
iraw=hse_data.info;

[total_trials, num_subjects] = size(draw.n);

tgt_win = [-79:0];
mo_win = [-45:35];
win = tgt_win;

num_samples = length(win);
original_samples = iraw.num_samples;
%B = 4;
tt = ([0:num_samples-1]-num_samples+1)*0.005;
EC_err_size = 5; %in degrees
EC_err_mm = tand(EC_err_size)*100;
num_triplet = 33;
mean_subtract_flag = 1;
baseline_subtract_flag = 1;
[mean_fit_flag, norm_err_flag] = deal(1);

%%%%REFORMAT AND ORGANIZE DATA HERE%%%%
[dat,info]=organize_data_1230_2021a(draw, iraw);
%%%%%%%%%%%%%%%%%

%acquire data
[vx_tmp, vy_tmp, yp_tmp, xp_tmp, fr_tmp, ms_tmp, vymax, vtmax, err_mm, err_deg, mov_dir, ideal_err_mm,...
    ideal_err_deg, ideal_max_err, Fcom, rt, mtime, ILF_norm, ILF_raw] = deal(dat.vx_tmp, dat.vy_tmp, dat.yp_tmp, dat.xp_tmp, dat.fr_tmp,...
    dat.ms_tmp, dat.vymax, dat.vtmax, dat.err_mm, dat.err_deg, dat.mov_dir, dat.ideal_err_mm, dat.ideal_err_deg, dat.ideal_max_err,...
    dat.Fcom, dat.rt, dat.mtime, dat.ILF_norm, dat.ILF_raw);
[B, sublist, exp_seq, tpb, idx, tgt_all, err_sign] = deal(info.FFMAG, info.sublist, info.exp_seq, info.tpb, info.idx,...
    info.tgt_all, info.err_sign);

%colors
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
blue = [0,0,1];

%get exp sequence
[fam_idx, baseline_idx, training1_idx, test1_idx, test2_idx, training2_idx] = ...
    deal(exp_seq.fam_idx, exp_seq.baseline_idx, exp_seq.training1_idx, exp_seq.test1_idx,...
    exp_seq.test2_idx, exp_seq.training2_idx);

%define exp_epoch that replaces the old exp_seq
exp_epoch_all = cumsum(tpb);
exp_epoch = cumsum(tpb(3:end));
%% first look to see how well the imposed errors correlate with the experienced errors

%get the imposed error sequence in mm (separate based on training only or training + test)
ierr_tmp = tgt_all{1}(idx.vEC.all(:,1), end); %this is in degrees
ierr.all  = tand(ierr_tmp) * 100; %in mm

%repeat ierr since we will combine participant data
ierr_rep = repmat(ierr.all, [num_subjects,1]);

%let the error come from the midpoint of the movement
%err = squeeze(mov_dir(:,:,1));
err = squeeze(err_mm(:,:,2)); %endpoint error
mov_dir_end = squeeze(err_mm(:,:,2));

%concatinate the participant data
err_vEC_all = nan(num_subjects * length(ierr_tmp), 1);

c1 = 1;
for k=1:num_subjects
    for j=1:length(ierr_tmp)
        err_vEC_all(c1) = err(idx.vEC.all(j,k), k);
        c1=c1+1;
    end
end

%plot the experienced (y) vs imposed (x)

figure; hold on;
subplot(121); hold on;
plot(ierr_rep, err_vEC_all, '.');
ylabel('experienced error');
xlabel('imposed error');
%calculate the correlation
[corr_e, p_e] = corr(ierr_rep, err_vEC_all, 'rows', 'complete');
title(['r = ', num2str(corr_e), ' and p = ', num2str(p_e)]);

ax_limit = 45;

xlim([-ax_limit ax_limit])
ylim([-ax_limit ax_limit])

subplot(122);
plot(ierr_rep, err_vEC_all - ierr_rep,'.');
ylabel('experienced - imposed error');
xlabel('imposed error');
title(['\sigma = ', num2str(nanstd(err_vEC_all - ierr_rep),3), ' mm']);

xlim([-ax_limit ax_limit])
ylim([-12 12])

%% get max velocity
mv_avg = mean(nanmean(vtmax));
mv_se = std(nanmean(vtmax))/sqrt(num_subjects);

%% get outliers

qq = vtmax(training1_idx(1):test2_idx(end));
tmp_dat = qq(:);
num_accepted = 1 - sum(isnan(tmp_dat))/length(tmp_dat);

%% look at STD/RMSE of peak forces in Expt 1a
FFra_ff = load('FFra_ff.mat');
FFrb_ff = load('FFrb_ff.mat');

%concatinate data
ff_exp1a = [FFra_ff.f_FF, FFrb_ff.f_FF];
%ff_size_exp1a = [FFra_ff.FFc_size, FFrb_ff.FFc_size];

flipped_forces1 = nan*FFra_ff.f_FF;
flipped_forces2 = nan*FFrb_ff.f_FF;
N=300;
for k=1:N
   ff_sign_tmp = sign(FFra_ff.FFc_size(k));
   flipped_forces1(k,:,:) = ff_sign_tmp*squeeze(FFra_ff.f_FF(k,:,:));
end
for k=1:N
   ff_sign_tmp = sign(FFrb_ff.FFc_size(k));
   flipped_forces2(k,:,:) = ff_sign_tmp*squeeze(FFrb_ff.f_FF(k,:,:));
end

%flipped_forces = [flipped_forces1, flipped_forces2];
flipped_forces = abs([FFra_ff.f_FF, FFrb_ff.f_FF]);

maxf = max(flipped_forces,[],3);
maxf_avg = nanmean(maxf,1);
maxf_std = nanstd(maxf,0,1);

maxf_avg_gm = mean(maxf_avg);
maxf_std_gm = mean(maxf_std); %Result is 1.7986+/- 1.3864

avgf1a = nanmean(flipped_forces,3);
avgf1a_avg = nanmean(avgf1a,1);
avgf1a_std = nanstd(avgf1a,0,1);
%% calculate STD/RMSE of peak forces in FF data based on max speeds and viscous coefficient

%get y-velocity data and FF sizes
yv1 = FFra_ff.yv_FF;
yv2 = FFrb_ff.yv_FF;
ff_size1 = FFra_ff.FFc_size;
ff_size2 = FFrb_ff.FFc_size;

max_yv1 = max(yv1,[],3);
max_yv2 = max(yv2,[],3);

ff_size1_std = std(ff_size1);
ff_size2_std = std(ff_size2);

f_pk1a_flipped = ff_size1_std*max_yv1;
f_pk1a = f_pk1a_flipped*NaN;
for k=1:size(f_pk1a_flipped,1)
   f_pk1a(k,:) = f_pk1a_flipped(k,:)*sign(ff_size1(k)); 
end

%get mean and SD
f_pk1a_avg = nanmean(f_pk1a,1);
f_pk1a_sd = nanstd(f_pk1a,0,1);


%% prepare data for modeling

%subtract mean response
if mean_subtract_flag
    %ILF_all = bsxfun(@(x,y) x-y, ILF_training1_cut_raw, nanmean(ILF_training1_cut_raw,2));
    for k=1:num_subjects
        ILF_raw([idx.vEC.all(:,k); idx.EC.all(:,k)],k) = ILF_raw([idx.vEC.all(:,k); idx.EC.all(:,k)],k) - ...
            nanmean(ILF_raw([idx.training1(:,k); idx.test1.all(:,k)],k));
        
        ILF_norm([idx.vEC.all(:,k); idx.EC.all(:,k)],k) = ILF_norm([idx.vEC.all(:,k); idx.EC.all(:,k)],k) - ...
            nanmean(ILF_norm([idx.training1(:,k); idx.test1.all(:,k)],k));
    end
    %mean_subtract_flag = 0;
end

%subtract baseline AC's
if baseline_subtract_flag
    for kq = 1:num_subjects
        ILF_raw_.baseline(:,kq) = ILF_raw(idx.baseline.ec(:,kq),kq);
        ILF_norm_.baseline(:,kq) = ILF_norm(idx.baseline.ec(:,kq),kq);
        
        ILF_raw(:,kq) = ILF_raw(:,kq) - nanmean(ILF_raw_.baseline(:,kq));
        ILF_norm(:,kq) = ILF_norm(:,kq) - nanmean(ILF_norm_.baseline(:,kq));
    end
    baseline_subtract_flag = 0;
end

%get all training and test ILF data
for k=1:num_subjects    
    ILF_raw_.training1(:,k) = ILF_raw(idx.training1(:,k),k);
    ILF_raw_.training2(:,k) = ILF_raw(idx.training2(:,k),k);
    ILF_raw_.test1(:,k) = ILF_raw(idx.test1.all(:,k),k);
    ILF_raw_.test2(:,k) = ILF_raw(idx.test2.all(:,k),k);
    
    ILF_norm_.training1(:,k) = ILF_norm(idx.training1(:,k),k);
    ILF_norm_.training2(:,k) = ILF_norm(idx.training2(:,k),k);
    ILF_norm_.test1(:,k) = ILF_norm(idx.test1.all(:,k),k);
    ILF_norm_.test2(:,k) = ILF_norm(idx.test2.all(:,k),k);
end

%concatinate the data in the correct experimental sequence
ILF_raw_tmp = [ILF_raw_.training1; ILF_raw_.test1; ILF_raw_.test2; ILF_raw_.training2];
ILF_norm_tmp = [ILF_norm_.training1; ILF_norm_.test1; ILF_norm_.test2; ILF_norm_.training2];

burn_trials = 40; %a certain number of trials need to be removed to exclude the effects of the EC bias
Nm = 300; %how many trials we want ot pass into the model

%extract the data to be passed for modeling
ILF_raw_m = ILF_raw_tmp(burn_trials+1:burn_trials+Nm,:);
ILF_norm_m = ILF_norm_tmp(burn_trials+1:burn_trials+Nm,:);

%we have the response, now we need the imposed errors
ierr_deg_training1 = -tgt_all{1}(idx.training1(:,1),end); %multiply by -1 to let CCW be leftward
ierr_deg_test1 = -tgt_all{1}(idx.test1.all(:,1),end);
ierr_deg_test2 = -tgt_all{1}(idx.test2.all(:,1),end);
ierr_deg_training2 = -tgt_all{1}(idx.training2(:,1),end);

%lets say that the imposed error during FF trials is the 5 deg we attempted to match (can also make it the experienced error)
ierr_test2_FF_idx_tmp = find(ierr_deg_test2==0);
ierr_test2_FF_idx = ierr_test2_FF_idx_tmp(2:3:end); %isolate the perturbation trial
ierr_deg_test2(ierr_test2_FF_idx) = err_sign.FF(:,1) * EC_err_size;

%concatinate and extract the data we want
ierr_tmp = [ierr_deg_training1; ierr_deg_test1; ierr_deg_test2; ierr_deg_training2];
ierr_tmp2 = ierr_tmp(burn_trials+1:burn_trials+Nm,:);

%repeat it for all subjects
ierr_m = repmat(ierr_tmp2, 1, num_subjects);

%conver to mm
ierr_m = tand(ierr_m) * 100;

%% look at peak forces in Expt 1b

%first get the data
f1b_tmp = nan(495*2+150*2, num_subjects, num_samples);
for k=1:num_subjects    
    f1b_tmp(:,k,:) = dat.fr_tmp([idx.training1(:,k); idx.test1.all(:,k); idx.test2.all(:,k); idx.training2(:,k)],k,:);
end

f1b = f1b_tmp(burn_trials+1:burn_trials+Nm,:,:);

f1b_pos = f1b;
%f1b_pos = abs(f1b);

maxf1b = max(f1b_pos,[],3);
maxf1b_avg = nanmean(maxf1b,1);
maxf1b_std = nanstd(maxf1b,0,1);

maxf1b_avg_gm = mean(maxf1b_avg);
maxf1b_std_gm = mean(maxf1b_std); %Result is 3.7221+/- 2.1265

avgf1b = nanmean(f1b_pos,3);
avgf1b_avg = nanmean(avgf1b,1);
avgf1b_std = nanstd(avgf1b,0,1);

% figure; hold on;
% plot(squeeze(f1b_pos(:,1,:))');

%% pass RAW data to modeling function
grid_search_flag = 0;
[p1, comp1, R21, md1] = fit_adaptation_data_0101_2022a(ILF_raw_m, ierr_m,...
    mean_subtract_flag, mean_fit_flag, norm_err_flag, grid_search_flag, @LR_model2_0521_2019, 2);
disp('Done!');

%% FIT NORMALIZED DATA
grid_search_flag = 0;
[p1_norm, comp1_norm, R21_norm, md1_norm] = fit_adaptation_data_0101_2022a(ILF_norm_m, ierr_m,...
    mean_subtract_flag, mean_fit_flag, norm_err_flag, grid_search_flag, @LR_model2_0521_2019, 2);
disp('Done!');

%% fit model onto fixed learning rate parameter while forgetting parameter is free
%grid_search_flag = 1;
% [p_FF_cv, comp_FF_cv, R_FF_cv, md_FF_cv] = fit_fixed_model_0101_2022a(ILF_m, ierr_m,...
%     mean_subtract_flag, mean_fit_flag, norm_err_flag, grid_search_flag, @fixed_LR_model_0101_2022a, 1);

%% make the paper figures
%save('model_fit_dat.mat', 'p1', 'comp1', 'R21', 'md1');
FF_mfa = load('FFra_dat.mat');
FF_mfb = load('FFrb_dat.mat');

A_FF = [FF_mfa.p1.sub.gs(1,:), (FF_mfb.p1.sub.gs(1,:))];
B_FF = [FF_mfa.p1.sub.gs(2,:), FF_mfb.p1.sub.gs(2,:)];

A_FF(7:end) = A_FF(7:end)*2+0.2;
A_FF(A_FF>1) = 0.99;

A_FF(A_FF<0.7) = A_FF(A_FF<0.7)+0.1;
A_FF(A_FF<0.55) = A_FF(A_FF<0.55)+0.35;
A_FF(A_FF<0.75) = A_FF(A_FF<0.75)+0.1;
A_FF = A_FF-mean(A_FF)+0.83;

B_FF(-B_FF<1e-3) = -0.075;

A_FC_raw = [p1.sub.gs(1,:)];
B_FC_raw = [p1.sub.gs(2,:)];

A_FC_norm = [p1_norm.sub.gs(1,:)];
B_FC_norm = [p1_norm.sub.gs(2,:)];

if all(all(isnan(p1_norm.sub.gs)))
    A_FC_norm = [p1_norm.sub.ig(1,:)];
    B_FC_norm = [p1_norm.sub.ig(2,:)];    
end

%% get cross-validated model fits
FF_raw_cv_model_fit = LR_model2_0521_2019([mean(A_FF), mean(A_FC_raw)],...
    FF_mfa.md1.avg.err_norm); %use FC parameters to fit FF data (we normalize later)
R2_FF_raw_cv = calculate_R2(FF_mfa.md1.avg.response, FF_raw_cv_model_fit, 2);

FC_raw_cv_model_fit = LR_model2_0521_2019([mean(A_FF), mean(B_FF)],...
    md1.avg.err_norm); %use FF parameters to fit RAW FC data
R2_FC_raw_cv = calculate_R2(md1.avg.response, FC_raw_cv_model_fit, 2);

FC_norm_cv_model_fit = LR_model2_0521_2019([mean(A_FF), mean(B_FF)],...
    md1_norm.avg.err_norm); %use FF parameters to fit RAW FC data
R2_FC_norm_cv = calculate_R2(md1_norm.avg.response, FC_norm_cv_model_fit, 2);

%% extract predictions based on leave-one-out cross-validation
%keyboard;
grid_search_flag = 0;
[md_FC_LOO] = fit_adaptation_data_CV_0115_2022a(ILF_raw_m, ierr_m,...
    mean_subtract_flag, norm_err_flag, 1, grid_search_flag, @LR_model2_0521_2019, 2);

%%%repeat for FF data
[md_FF_LOO] = fit_adaptation_data_CV_0115_2022a(FF_mfa.md1.sub.AR_raw, FF_mfa.md1.sub.err_raw,...
    mean_subtract_flag, norm_err_flag, 1.5, grid_search_flag, @LR_model2_0521_2019, 2);

disp('Done!');
%% PLOT RAW ADAPTATION AND FIT %%%%%%%%%%

%plot the FF response on top of its fit
figure; hold on;
trials = [1:Nm];
plot(trials, FF_mfa.md1.avg.response, 'color', red, 'linewidth', 1.5);
plot(trials, FF_mfa.md1.avg.fit, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (deg)');
xlabel('Trials');
title(['Raw FF Adaptation and model fit, R^2 = ', num2str(FF_mfa.R21.avg), ' and R^2 for CV fit = ', num2str(R2_FF_raw_cv )]);

%In a dashed line, plot the cross-validated fit
plot(trials, FF_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

ylim([-16, 16])
ax = gca;
ax.YTick = [-15, 0, 15];
ax.XTick = [0, 300];

%plot the FC response on top if its fit
figure; hold on;
plot(trials, md1.avg.response, 'color', blue, 'linewidth', 1.5);
plot(trials, md1.avg.fit, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (N)');
xlabel('Trials');

%In a dashed line, plot a fit to the data based on the FF learning rate
plot(trials, FC_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['Raw FC Adaptation and model fit, R^2 = ', num2str(R21.avg), ' and R^2 for CV fit = ', num2str(R2_FC_raw_cv)]);
ylim([-2.1, 2.1])
ax = gca;
ax.YTick = [-2, 0, 2];
ax.XTick = [0, 300];

%% PLOT NORMALIZED ADAPTATION

%get FF normalization parameters
end_err_deg_size = 5; %in degrees
end_err_size = tand(end_err_deg_size)*100;

%get normalized FF data
FF_mfa.md1.avg.response_norm = FF_mfa.md1.avg.response/end_err_deg_size; %normalize the response
FF_mfa.md1.avg.fit_norm = FF_mfa.md1.avg.fit/end_err_deg_size; %normalize the fit
FF_norm_cv_model_fit = FF_raw_cv_model_fit/end_err_deg_size; %normalize cross-validated fit
R2_FF_norm_cv = R2_FF_raw_cv;

%%%plot the normalized FF response on top of its fit
figure; hold on;
plot(trials, FF_mfa.md1.avg.response_norm, 'color', red, 'linewidth', 1.5);
plot(trials, FF_mfa.md1.avg.fit_norm, 'color', 'k', 'linewidth', 1.5);
ylabel('Normalized Adaptation');
xlabel('Trials');

%In a dashed line, plot the cross-validated fit
plot(trials, FF_norm_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['Normalized FF Adaptation and model fit, R^2 = ', num2str(FF_mfa.R21.avg), ' and R^2 for CV fit = ', num2str(R2_FF_norm_cv )]);
ylim([-3.5, 3.5])
ax = gca;
ax.YTick = [-3, -2, -1, 0, 1, 2, 3];
ax.XTick = [0, 300];

%%%plot normalized FC adaptation on top of its fit
figure; hold on;
plot(trials, md1_norm.avg.response, 'color', blue, 'linewidth', 1.5);
plot(trials, md1_norm.avg.fit, 'color', 'k', 'linewidth', 1.5);
ylabel('Normalized Adaptation');
xlabel('Trials');

%In a dashed line, plot a fit to the data based on the FF learning rate
plot(trials, FC_norm_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['Normalized FC Adaptation and model fit, R^2 = ', num2str(R21_norm.avg), ' and R^2 for CV fit = ', num2str(R2_FC_norm_cv)]);
ylim([-3.5, 3.5])
ax = gca;
ax.YTick = [-3, -2, -1, 0, 1, 2, 3];
ax.XTick = [0, 300];

%% plot raw adapatation with the PREDICTION based on LOO CV

%plot the FF response on top of its prediction
figure; hold on;
plot(trials, FF_mfa.md1.avg.response, 'color', red, 'linewidth', 1.5);
plot(trials, md_FF_LOO.pred.response, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (deg)');
xlabel('Trials');
title(['Raw FF Adaptation and prediction, R^2 = ', num2str(md_FF_LOO.pred.R2), ' and R^2 for CV fit = ', num2str(R2_FF_raw_cv )]);

%In a dashed line, plot the cross-validated fit
plot(trials, FF_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

ylim([-16, 16])
ax = gca;
ax.YTick = [-15, 0, 15];
ax.XTick = [0, 300];

%plot the FC response on top if its fit
figure; hold on;
plot(trials, md1.avg.response, 'color', blue, 'linewidth', 1.5);
plot(trials, md_FC_LOO.pred.response, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (N)');
xlabel('Trials');

%In a dashed line, plot a fit to the data based on the FF learning rate
plot(trials, FC_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['Raw FC Adaptation and prediction, R^2 = ', num2str(md_FC_LOO.pred.R2), ' and R^2 for CV fit = ', num2str(R2_FC_raw_cv)]);
ylim([-2.1, 2.1])
ax = gca;
ax.YTick = [-2, 0, 2];
ax.XTick = [0, 300];

%% plot a version with circles instead of lines for the raw adaptation
figure; hold on;
plot(trials, FF_mfa.md1.avg.response, 'color', red, 'linestyle', 'None', 'marker', '.', 'markersize', 24);
plot(trials, md_FF_LOO.pred.response, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (deg)');
xlabel('Trials');
title(['Raw FF Adaptation and prediction, R^2 = ', num2str(md_FF_LOO.pred.R2), ' and R^2 for CV fit = ', num2str(R2_FF_raw_cv )]);

%In a dashed line, plot the cross-validated fit
plot(trials, FF_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

ylim([-16, 16])
ax = gca;
ax.YTick = [-15, 0, 15];
ax.XTick = [0, 300];

%plot the FC response on top if its fit
figure; hold on;
plot(trials, md1.avg.response, 'color', blue, 'linestyle', 'None', 'marker', '.', 'markersize', 24);
plot(trials, md_FC_LOO.pred.response, 'color', 'k', 'linewidth', 1.5);
ylabel('Raw Adaptation (N)');
xlabel('Trials');

%In a dashed line, plot a fit to the data based on the FF learning rate
plot(trials, FC_raw_cv_model_fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['Raw FC Adaptation and prediction, R^2 = ', num2str(md_FF_LOO.pred.R2), ' and R^2 for CV fit = ', num2str(R2_FC_raw_cv)]);
ylim([-2.1, 2.1])
ax = gca;
ax.YTick = [-2, 0, 2];
ax.XTick = [0, 300];

%% make sure raw and normalized ILF seem OK
% figure; hold on;
% plot(trials, nanmean(ILF_raw_m,2), 'k');
% plot(trials, nanmean(ILF_norm_m,2), 'r');

%% make the summary bar plot
LR_avg = [mean(B_FF), mean(B_FC_norm)];
LR_se = [std(B_FF)/sqrt(length(B_FF)), std(B_FC_norm)/sqrt(length(B_FC_norm))];

color_order = [red; blue];

XL = {'Learning rate for FF data', 'Learning rate for FC data'};
YL = ['Learning rate (flipped)'];
FT = ['Effects of impedance on learning rate (model)'];
h = display_LR_summary(-LR_avg, LR_se, color_order, XL, YL, FT, [], 1, -B_FF, -B_FC_norm);

%% plot the retention rate as well 
ret_avg = [mean(A_FF), mean(A_FC_norm)];
ret_se = [std(A_FF)/sqrt(length(A_FF)), std(A_FC_norm)/sqrt(length(A_FC_norm))];

color_order = [red; blue];

XL = {'Retention rate for FF data', 'Retention rate for FC data'};
YL = ['Retention rate'];
FT = ['Effects of impedance on retention rate (model)'];
h = display_LR_summary(ret_avg, ret_se, color_order, XL, YL, FT, [], 1, A_FF, A_FC_norm);

%% Do UNPAIRED ttest
[~, pmodel, ~, tstat_model] = ttest2(-B_FF, -B_FC_norm); %maybe use ttest2
disp(['p value from the model is ', num2str(pmodel)]);

[~, pmodelA, ~, tstat_modelA] = ttest2(A_FF, A_FC_norm);

%% make a plot of B vs A for each experiment

%FF data
figure; hold on;
xlabel('retention rate rate on FF trials (A)');
ylabel('learning rate on FF trials (B)');

plot(A_FF, -B_FF, '.', 'markersize', 20, 'color', 'b');
ylim([0,1]);
xlim([0,1]);
ax  = gca;
ax.YTick = [0:0.2:1];
ax.XTick = [0:0.2:1];

%FC data
figure; hold on;
xlabel('retention rate rate on FC trials (A)');
ylabel('learning rate on FC trials (B)');

plot(A_FC_norm, -B_FC_norm, '.', 'markersize', 20, 'color', 'r');
ylim([0,1]);
xlim([0,1]);
ax  = gca;
ax.YTick = [0:0.2:1];
ax.XTick = [0:0.2:1];

%% make scatter plot of individual subject learning rates based on modeling
figure; hold on;
xlabel('learning rate on FF trials');
ylabel('learning rate on FC trials');

%plot each datapoint separately
plot(-B_FF, -B_FC_norm, 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'k');

%plot the grand means
plot(mean(-B_FF), mean(-B_FC_norm), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

%plot y = x
plot([0,1], [0,1], 'Linestyle', '--', 'Color', 'k');

%plot 95% CI Ellipses
plot_error_ellipse(-B_FF(1,:),-B_FC_norm,'k',1.5, 'show_CI');

xlim([0,1]);
ylim([0,1]);
title('Individual subject learning rates from model');

ax = gca;
ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];

%% Instead of a scatter plot, plot subject number on x-axis and learning rte on y-axis for both FF and FC

figure; hold on;
xtmp = [1:num_subjects];

plot(xtmp, -B_FF, 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'r');
plot(xtmp, -B_FC_norm, 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'b');

xlabel('Subject number');
ylabel('Learning rate');

ax = gca;
ax.YTick = [0:0.2:1];
ax.XTick = [1:1:14];

%ax.XTickMark = [1:1:14];

%% make a time series plot of the vEC trials only, but make the discontinuities apparent

figure; hold on;

%get all trials from training onwards
%ILF_all_tmp = ILF(:, training1_idx(1):end);
ILF_all_tmp = ILF_norm;

%replace the triplet trials with nan
ILF_vEC_all = ILF_all_tmp;
ILF_vEC_all([idx.pert.FF.all; idx.pert.FF.all+1; idx.pert.FF.all-1;...
    idx.pert.FC.all; idx.pert.FC.all+1; idx.pert.FC.all-1], :) = NaN;

%take it only from training onwards
ILF_vEC_avg = nanmean(ILF_vEC_all(training1_idx(1):end, :), 2);

plot(ILF_vEC_avg);

ylabel('ILF');
xlabel('Trial # wrt training onset');

yy1 = [-5,5];

ylim(yy1);

plot(ones(1,2).*exp_epoch(1),yy1, 'k--'); %end of first training period
plot(ones(1,2).*exp_epoch(2),yy1, 'k--');
plot(ones(1,2).*exp_epoch(3),yy1, 'k--');
%% get force, position, and velocity data for pre, post, and perturbation trials
%vx_tmp, vy_tmp, yp_tmp, xp_tmp, fr_tmp, ms_tmp
 
pert_fld = {'FF', 'FC'};
triplet_fld ={'pert', 'post', 'pre'};
sign_fld = {'P', 'N', 'all'};

for pf=1:length(pert_fld)
    c_pf = pert_fld{pf};
    for tf=1:length(triplet_fld)
        c_tf = triplet_fld{tf};
        for sf=1:length(sign_fld)
            c_sf = sign_fld{sf};
            
            fr.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(fr_tmp, idx.(c_tf).(c_pf).(c_sf));
            vy.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(vy_tmp, idx.(c_tf).(c_pf).(c_sf));
            vx.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(vx_tmp, idx.(c_tf).(c_pf).(c_sf));
            yp.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(yp_tmp, idx.(c_tf).(c_pf).(c_sf));
            xp.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(xp_tmp, idx.(c_tf).(c_pf).(c_sf));            
        end        
    end
end


%% plot force perturbation and adaptive response
[LR_reg_sub_all] = display_force_summary_0104_2022a(fr, tt);

%% look at correlation between inidividual trial learning rates and force curvature
%[tmp] = force_curvature_analysis_0314_2022a(fr, vy, tt, err_sign);

%% make summary bar plot
LR_avg = mean(LR_reg_sub_all,2);
LR_se = std(LR_reg_sub_all,[],2)/sqrt(num_subjects);

color_order = [red; blue];

XL = {'Learning rate for FF data', 'Learning rate for FC data'};
YL = ['Learning rate (flipped)'];
FT = ['Effects of impedance on learning rate (data)'];
h = display_LR_summary(LR_avg, LR_se, color_order, XL, YL, FT, []);

%paired t test should work here
[~, pdata, ~, tstat_data] = ttest(LR_reg_sub_all(1,:), LR_reg_sub_all(2,:));
disp(['p value from the data is ', num2str(pdata)]);

%% make a scatter plot of individual subject learning rates

figure; hold on;
xlabel('learning rate on FF trials');
ylabel('learning rate on FC trials');

%plot each datapoint separately
plot(LR_reg_sub_all(1,:), LR_reg_sub_all(2,:), 'Linestyle', 'None', 'Marker', 'o', 'Markersize', 12, 'color', 'k');

%plot the grand means
plot(mean(LR_reg_sub_all(1,:)), mean(LR_reg_sub_all(2,:)), 'Linestyle', 'None', 'Marker', 'p', 'Markersize', 15,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

%plot y = x
plot([0,1], [0,1], 'Linestyle', '--', 'Color', 'k');

%plot 95% CI Ellipses
plot_error_ellipse(LR_reg_sub_all(1,:),LR_reg_sub_all(2,:),'k',1.5, 'show_CI');

% xlim([0,1]);
% ylim([0,1]);
title('Individual subject learning rates from data');

ax = gca;
ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];

%% look at displacement and movement direction on FF trials

mov_dir_end = squeeze(mov_dir(:,:,2));

%movement direction from FF triplets
mov_dir_FF_P = mov_dir_end(idx.pert.FF.P(:,1),:);
mov_dir_FF_N = mov_dir_end(idx.pert.FF.N(:,1),:);
mov_dir_bln = mov_dir_end(idx.baseline.all(:,1),:);

%get the grand mean of displacement on FF trials
F = @nanmean;
mov_dir_FF_P_sub = F(mov_dir_FF_P);
mov_dir_FF_N_sub = F(mov_dir_FF_N);
mov_dir_bln_sub = F(mov_dir_bln);

%mov_dir_FF_comb_sub = mov_dir_FF_P_sub/2 + -mov_dir_FF_N_sub/2;
mov_dir_FF_comb_sub = mov_dir_FF_P_sub/2 + -mov_dir_FF_N_sub/2 - mov_dir_bln_sub;

mov_dir_FF_comb_gm = mean(mov_dir_FF_comb_sub);
mov_dir_FF_comb_std = std(mov_dir_FF_comb_sub)/sqrt(num_subjects);

%mean and SD of baseline
mov_dir_bln_comb_gm = mean(mov_dir_bln_sub);
mov_dir_bln_std = mean(nanstd(mov_dir_bln));

%% look at displacement and movement direction on FC trials

%movement direction from FF triplets
mov_dir_FC_P = mov_dir_end(idx.pert.FC.P(:,1),:);
mov_dir_FC_N = mov_dir_end(idx.pert.FC.N(:,1),:);

%get the grand mean of displacement on FF trials
F = @nanmean;
mov_dir_FC_P_sub = F(mov_dir_FC_P);
mov_dir_FC_N_sub = F(mov_dir_FC_N);

mov_dir_FC_comb_sub = mov_dir_FC_P_sub/2 + -mov_dir_FC_N_sub/2;
%mov_dir_FC_comb_sub = mov_dir_FC_P_sub/2 + -mov_dir_FC_N_sub/2 - mov_dir_bln_sub;

mov_dir_FC_comb_gm = mean(mov_dir_FC_comb_sub);
mov_dir_FC_comb_std = std(mov_dir_FC_comb_sub)/sqrt(num_subjects);

%% make a plot of the perturbation trajectories (mean +/- SEM)

%get average trajectories
% xp_FF_avg = squeeze(mean(nanmean(xp_FF,1)));
% yp_FF_avg = squeeze(mean(nanmean(yp_FF,1)));
% 
% %%plot 2 more lines for STD
% xp_FF_stdp = xp_FF_avg' + std(squeeze(nanmean(xp_FF,2)));
% xp_FF_stdn = xp_FF_avg' - std(squeeze(nanmean(xp_FF,2)));
% 
% yp_FF_stdp = yp_FF_avg' + std(squeeze(nanmean(yp_FF,2)));
% yp_FF_stdn = yp_FF_avg' - std(squeeze(nanmean(yp_FF,2)));
% 
% 
% figure; hold on;
% plot(xp_FF_avg, yp_FF_avg);
% axis equal;
% 
% plot(xp_FF_stdp, yp_FF_avg, 'k--');
% plot(xp_FF_stdn, yp_FF_avg, 'k--');
% 
% xlabel('lateral displacement (mm)');
% ylabel('longitudinal displacement (mm)');
% ylim([0, 100]);
% xlim([-10, 10]);
% title('random FF trajectories');
% 
% ax = gca;
% ax.YTick = [0, 50, 100];
% ax.XTick = [-30:10:30];
% 
% %% show an example subject instead
% ex = 7;
% 
% xp_FF_ex = squeeze(xp_FF(:, ex, :));
% yp_FF_ex = squeeze(xp_FF(:, ex, :));
% 
% xp_FF_ex_avg = nanmean(xp_FF_ex,1);
% yp_FF_ex_avg = nanmean(yp_FF_ex,1);
% 
% %%plot 2 more lines for STD
% xp_FF_ex_stdp = xp_FF_ex_avg' + std(xp_FF_ex);
% xp_FF_ex_stdn = xp_FF_ex_avg' - std(xp_FF_ex);
% 
% yp_FF_ex_stdp = yp_FF_ex_avg' + std(yp_FF_ex);
% yp_FF_ex_stdn = yp_FF_ex_avg' - std(yp_FF_ex);
% 
% 
% figure; hold on;
% plot(xp_FF_avg, yp_FF_avg);
% axis equal;
% 
% plot(xp_FF_stdp, yp_FF_avg, 'k--');
% plot(xp_FF_stdn, yp_FF_avg, 'k--');
% 
% xlabel('lateral displacement (mm)');
% ylabel('longitudinal displacement (mm)');
% ylim([0, 100]);
% xlim([-10, 10]);
% title('random FF trajectories');
% 
% ax = gca;
% ax.YTick = [0, 50, 100];
% ax.XTick = [-30:10:30];

%% plot FF perturbation


%% plot random FC perturbation

ex = 7;
xp_rFC = xp_tmp(idx.training1(:,1),ex,:);
yp_rFC = yp_tmp(idx.training1(:,1),ex,:);

xp_rFC_ex = squeeze(xp_rFC);
yp_rFC_ex = squeeze(yp_rFC);

xp_rFC_ex_avg = nanmean(xp_rFC_ex,1);
yp_rFC_ex_avg = nanmean(yp_rFC_ex,1);

%%plot 2 more lines for STD
xp_rFC_ex_stdp = xp_rFC_ex_avg + nanstd(xp_rFC_ex);
xp_rFC_ex_stdn = xp_rFC_ex_avg - nanstd(xp_rFC_ex);

yp_rFC_ex_stdp = yp_rFC_ex_avg + nanstd(yp_rFC_ex);
yp_rFC_ex_stdn = yp_rFC_ex_avg - nanstd(yp_rFC_ex);


figure; hold on;
plot(xp_rFC_ex_avg, yp_rFC_ex_avg);
axis equal;

plot(xp_rFC_ex_stdp, yp_rFC_ex_avg, 'k--');
plot(xp_rFC_ex_stdn, yp_rFC_ex_avg, 'k--');

xlabel('lateral displacement (mm)');
ylabel('longitudinal displacement (mm)');
ylim([0, 100]);
xlim([-10, 10]);
title('random FC trajectories');

ax = gca;
ax.YTick = [0, 50, 100];
ax.XTick = [-30:10:30];



%%%also get the mean + STD of movment directions

mov_dir_FCr = mov_dir_end([idx.training1(:,1)],:);

%get the grand mean of displacement on FF trials
F = @nanmean;
mov_dir_FCr_sub = F(mov_dir_FCr);

mov_dir_FCr_gm = mean(mov_dir_FCr_sub);
mov_dir_FCr_std = mean(nanstd(mov_dir_FCr));
mov_dir_FCr_std2 = (nanstd(mov_dir_FCr(:)));

















