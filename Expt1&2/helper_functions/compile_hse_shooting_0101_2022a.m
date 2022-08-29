clear all;
close all;
home;

%cd('C:\Users\laith\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\mat_files\ALL_DATA'); %laptop
cd('C:\Users\ryanm\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\mat_files\ALL_DATA'); %Ryan's computer!

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
hse_data=load(strcat(Cdr,'\',hse_filename));
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
    ideal_err_deg, ideal_max_err, Fcom, rt, mtime, ILF] = deal(dat.vx_tmp, dat.vy_tmp, dat.yp_tmp, dat.xp_tmp, dat.fr_tmp,...
    dat.ms_tmp, dat.vymax, dat.vtmax, dat.err_mm, dat.err_deg, dat.mov_dir, dat.ideal_err_mm, dat.ideal_err_deg, dat.ideal_max_err,...
    dat.Fcom, dat.rt, dat.mtime, dat.ILF);
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

%% prepare data for modeling

%subtract mean response
if mean_subtract_flag
    %ILF_all = bsxfun(@(x,y) x-y, ILF_training1_cut_raw, nanmean(ILF_training1_cut_raw,2));
    for k=1:num_subjects
        ILF([idx.vEC.all(:,k); idx.EC.all(:,k)],k) = ILF([idx.vEC.all(:,k); idx.EC.all(:,k)],k) - ...
            nanmean(ILF([idx.training1(:,k); idx.test1.all(:,k)],k));        
    end
    %mean_subtract_flag = 0;
end

%subtract baseline AC's
if baseline_subtract_flag
    for kq = 1:num_subjects
        ILF_.baseline(:,kq) = ILF(idx.baseline.ec(:,kq),kq);
        ILF(:,kq) = ILF(:,kq) - nanmean(ILF_.baseline(:,kq));
    end
    baseline_subtract_flag = 0;
end

%get all training and test ILF data
for k=1:num_subjects    
    ILF_.training1(:,k) = ILF(idx.training1(:,k),k);
    ILF_.training2(:,k) = ILF(idx.training2(:,k),k);
    ILF_.test1(:,k) = ILF(idx.test1.all(:,k),k);
    ILF_.test2(:,k) = ILF(idx.test2.all(:,k),k);
end

%concatinate the data in the correct experimental sequence
ILF_tmp = [ILF_.training1; ILF_.test1; ILF_.test2; ILF_.training2];

burn_trials = 40; %a certain number of trials need to be removed to exclude the effects of the EC bias
Nm = 300; %how many trials we want ot pass into the model

%extract the data to be passed for modeling
ILF_m = ILF_tmp(burn_trials+1:burn_trials+Nm,:);

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

%% pass vEC data to modeling function
grid_search_flag = 1;
[p1, comp1, R21, md1] = fit_adaptation_data_0101_2022a(ILF_m, ierr_m,...
    mean_subtract_flag, mean_fit_flag, norm_err_flag, grid_search_flag, @LR_model2_0521_2019, 2);
disp('Done!');

%% do a fit where the FF learning is used, but the forgetting parameter is free
grid_search_flag = 1;
[p2, comp2, R22, md2] = fit_fixed_model_0101_2022a(ILF_m, ierr_m,...
    mean_subtract_flag, mean_fit_flag, norm_err_flag, grid_search_flag, @fixed_LR_model_0101_2022a, 1);
disp('Done!');

%% make the paper figures
%save('model_fit_dat.mat', 'p1', 'comp1', 'R21', 'md1');
FF_mfa = load('model_fit_dat.mat');
FF_mfb = load('FFr_dat.mat');

A_FF = [FF_mfa.p1.sub.gs(1,:), FF_mfb.p1.sub.gs(1,:)];
B_FF = [FF_mfa.p1.sub.gs(2,:), FF_mfb.p1.sub.gs(2,:)];

A_FC = [p1.sub.gs(1,:), p1.sub.gs(1,:)];
B_FC = [p1.sub.gs(2,:)];

A_FC2 = [p2.sub.gs(:)];

%plot the FF response on top of its fit
figure; hold on;
trials = [1:Nm];
plot(trials, FF_mfa.md1.avg.response, 'color', orange, 'linewidth', 1.5);
plot(trials, FF_mfa.md1.avg.fit, 'color', 'k', 'linewidth', 1.5);
ylabel('Normalized Adaptation');
xlabel('Trials');
title(['FF Adaptation and model fit, R^2 = ', num2str(FF_mfa.R21.avg)]);
ylim([-2.5, 2.5])

%plot the FC response on top if its fit
figure; hold on;
plot(trials, md1.avg.response, 'color', orange, 'linewidth', 1.5);
plot(trials, md1.avg.fit, 'color', 'k', 'linewidth', 1.5);
ylabel('Normalized Adaptation');
xlabel('Trials');

%In a dashed line, plot a fit to the data based on the FF learning rate
plot(trials, md2.avg.fit, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');

title(['FC Adaptation and model fit, R^2 = ', num2str(R21.avg), ' and R^2 for B_F_F fit = ', num2str(R22.avg)]);
ylim([-2.5, 2.5])

%make the bar plot
LR_avg = [mean(B_FF), mean(B_FC)];
LR_se = [std(B_FF)/sqrt(length(B_FF)), std(B_FC)/sqrt(length(B_FC))];

color_order = [orange; cyan];

XL = {'Learning rate for FF data', 'Learning rate for FC data'};
YL = ['Learning rate (flipped)'];
FT = ['Effects of impedance on learning rate (model)'];
h = display_LR_summary(-LR_avg, LR_se, color_order, XL, YL, FT, []);

%Do UNPAIRED ttest
[~, pmodel] = ttest2(-B_FF, -B_FC);
disp(['p value from the model is ', num2str(pmodel)]);

%% make a time series plot of the vEC trials only, but make the discontinuities apparent

figure; hold on;

%get all trials from training onwards
%ILF_all_tmp = ILF(:, training1_idx(1):end);
ILF_all_tmp = ILF;

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
            yp.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(yp_tmp, idx.(c_tf).(c_pf).(c_sf));
            xp.(c_tf).(c_pf).(c_sf) = get_sub_data_0104_2022a(xp_tmp, idx.(c_tf).(c_pf).(c_sf));            
        end        
    end
end


%% plot force perturbation and adaptive response
[LR_reg_sub_all] = display_force_summary_0104_2022a(fr, tt);

%% make summary bar plot
LR_avg = mean(LR_reg_sub_all,2);
LR_se = std(LR_reg_sub_all,[],2)/sqrt(num_subjects);

color_order = [orange; cyan];

XL = {'Learning rate for FF data', 'Learning rate for FC data'};
YL = ['Learning rate (flipped)'];
FT = ['Effects of impedance on learning rate (data)'];
h = display_LR_summary(LR_avg, LR_se, color_order, XL, YL, FT, []);

%paired t test should work here
[~, pdata] = ttest(LR_reg_sub_all(1,:), LR_reg_sub_all(2,:));
disp(['p value from the data is ', num2str(pdata)]);

%% look at displacement and movement direction on FF trials

mov_dir_end = squeeze(mov_dir(:,:,2));

%movement direction from FF triplets
mov_dir_FF_P = mov_dir_end(idx.pert.FF.P(:,1),:);
mov_dir_FF_N = mov_dir_end(idx.pert.FF.N(:,1),:);

%get the grand mean of displacement on FF trials
mov_dir_FF_P_sub = nanmedian(mov_dir_FF_P);
mov_dir_FF_N_sub = nanmedian(mov_dir_FF_N);

mov_dir_FF_comb_sub = mov_dir_FF_P_sub/2 + mov_dir_FF_N_sub/2;

mov_dir_FF_comb_gm = mean(mov_dir_FF_comb_sub);
mov_dir_FF_comb_std = std(mov_dir_FF_comb_sub);




