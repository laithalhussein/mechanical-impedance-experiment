function [LR_reg_sub_all, ar_sub_all, ideal_avg_sub_all, ideal_median_sub_all,...
    ar_gm, ideal_gm] = display_force_summary_10_13_2021a(f, pert_fld, env_fld, sub_id, tt, show_median, show_error, show_PN, flip_sign)

%displays positive and negative force profiles for each
%participant, with the ideal, pre, post, and AR profiles

%Note that we display the difference between each case in a separate
%function
%also returns the adaptive response by regression

%thin blue = pre +, thin red = post +
%thick blue = pre -, thick red = pre +

force_UB = 6;
force_LB = -6;
grey = [0.5,0.5,0.5];
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
col_gain = 0.7;

[~, num_subjects, num_samples] = size(f.pert.(pert_fld).P);

%we will save the learning rate measurements, and the subject averaged
%ARs and ideal force data
LR_reg_sub_all = nan(num_subjects, 1);

ar_sub_all = nan(num_samples, num_subjects);
ideal_avg_sub_all = nan(num_samples, num_subjects);
ideal_median_sub_all = nan(num_samples, num_subjects);
pre_sub_all = nan(num_samples, num_subjects);
post_sub_all = nan(num_samples, num_subjects);

%%%make sure to have an extra subplot for the grand mean
if mod(num_subjects,2)==0,
    num_plots = ceil(num_subjects/2) + 1;
else
    num_plots = ceil(num_subjects/2);
end

%%%get the data
cf_pert_P = f.pert.(pert_fld).P;
cf_pert_N = f.pert.(pert_fld).N;

cf_pre_P = f.pre.(pert_fld).P;
cf_pre_N = f.pre.(pert_fld).N;

cf_post_P = f.post.(pert_fld).P;
cf_post_N = f.post.(pert_fld).N;

%calculate the grand mean of post, pre, ideal and AR trials
post_P_gm_tmp = nanmean(squeeze(nanmean(cf_post_P,1)),1);
post_N_gm_tmp = nanmean(squeeze(nanmean(cf_post_N,1)),1);
post_gm = post_P_gm_tmp/2 + -post_N_gm_tmp/2;

pre_P_gm_tmp = nanmean(squeeze(nanmean(cf_pre_P,1)),1);
pre_N_gm_tmp = nanmean(squeeze(nanmean(cf_pre_N,1)),1);
pre_gm = pre_P_gm_tmp/2 + -pre_N_gm_tmp/2;

ideal_P_gm_tmp = nanmean(squeeze(nanmean(cf_pert_P,1)),1);
ideal_N_gm_tmp = nanmean(squeeze(nanmean(cf_pert_N,2)),1);
ideal_gm = ideal_P_gm_tmp/2 + -ideal_N_gm_tmp/2;

%calculate the grand mean first so that we can put it on the individual subject plots
cf_AR_P_gm = (mean(squeeze(nanmean(cf_post_P,1)),1) - mean(squeeze(nanmean(cf_pre_P,1)),1));
cf_AR_N_gm = (mean(squeeze(nanmean(cf_post_N,1)),1) - mean(squeeze(nanmean(cf_pre_N,1)),1));
ar_gm = cf_AR_P_gm/2 + -cf_AR_N_gm/2;

%calculate standard errors for population averages
post_P_gm_err = nanstd(squeeze(nanmean(cf_post_P,1)),0,1)/sqrt(num_subjects);
post_N_gm_err = nanstd(squeeze(nanmean(cf_post_N,1)),0,1)/sqrt(num_subjects);
post_gm_tmp = squeeze(nanmean(cf_post_P,1))/2 + -squeeze(nanmean(cf_post_N,1))/2;
post_gm_err= nanstd(post_gm_tmp,0,1)/sqrt(num_subjects);

pre_P_gm_err = nanstd(squeeze(nanmean(cf_pre_P,1)),0,1)/sqrt(num_subjects);
pre_N_gm_err = nanstd(squeeze(nanmean(cf_pre_N,1)),0,1)/sqrt(num_subjects);
pre_gm_tmp = squeeze(nanmean(cf_pre_P,1))/2 + -squeeze(nanmean(cf_pre_N,1))/2;
pre_gm_err= nanstd(pre_gm_tmp,0,1)/sqrt(num_subjects);

ideal_P_gm_err = nanstd(squeeze(nanmean(cf_pert_P,1)),0,1)/sqrt(num_subjects);
ideal_N_gm_err = nanstd(squeeze(nanmean(cf_pert_N,1)),0,1)/sqrt(num_subjects);
ideal_gm_tmp = squeeze(nanmean(cf_pert_P,1))/2 + -squeeze(nanmean(cf_pert_N,1))/2;
ideal_gm_err= nanstd(ideal_gm_tmp,0,1)/sqrt(num_subjects);

ar_P_gm_tmp = squeeze(nanmean(cf_post_P,1)) - squeeze(nanmean(cf_pre_P,1));
ar_N_gm_tmp = squeeze(nanmean(cf_post_N,1)) - squeeze(nanmean(cf_pre_N,1));
ar_comb_gm_tmp = ar_P_gm_tmp/2 + -ar_N_gm_tmp/2;
ar_comb_gm_err = nanstd(ar_comb_gm_tmp,0,1)/sqrt(num_subjects);

%% repeat the above, but this time also calculate the median
post_P_med_tmp = nanmean(squeeze(nanmedian(cf_post_P,1)),1);
post_N_med_tmp = nanmean(squeeze(nanmedian(cf_post_N,1)),1);
post_med_gm = post_P_med_tmp/2 + -post_N_med_tmp/2;

pre_P_med_tmp = nanmean(squeeze(nanmedian(cf_pre_P,1)),1);
pre_N_med_tmp = nanmean(squeeze(nanmedian(cf_pre_N,1)),1);
pre_med_gm = pre_P_med_tmp/2 + -pre_N_med_tmp/2;

ideal_P_med_tmp = nanmean(squeeze(nanmedian(cf_pert_P,1)),1);
ideal_N_med_tmp = nanmean(squeeze(nanmedian(cf_pert_N,2)),1);
ideal_med_gm = ideal_P_med_tmp/2 + -ideal_N_med_tmp/2;

%calculate the grand mean first so that we can put it on the individual subject plots
cf_AR_P_med_gm = (mean(squeeze(nanmedian(cf_post_P,1)),1) - mean(squeeze(nanmedian(cf_pre_P,1)),1));
cf_AR_N_med_gm = (mean(squeeze(nanmedian(cf_post_N,1)),1) - mean(squeeze(nanmedian(cf_pre_N,1)),1));
ar_med_gm = cf_AR_P_med_gm/2 + -cf_AR_N_med_gm/2;


%% loop through subjects for plotting
%%%begin with plotting
figure; cc = 1;
for isub=1:num_subjects
    
    %determine the adaptive response and ideal force based on means
    cf_pert_P_sub_avg = nanmean(squeeze(cf_pert_P(:,isub,:)),1);
    cf_pert_N_sub_avg = nanmean(squeeze(cf_pert_N(:,isub,:)),1);
    cf_pert_comb_sub_avg = cf_pert_P_sub_avg/2 + -cf_pert_N_sub_avg/2;
    
    cf_post_P_sub_avg = nanmean(squeeze(cf_post_P(:,isub,:)),1);
    cf_post_N_sub_avg = nanmean(squeeze(cf_post_N(:,isub,:)),1);
    cf_post_comb_sub_avg = cf_post_P_sub_avg/2 + -cf_post_N_sub_avg/2;
    
    cf_pre_P_sub_avg = nanmean(squeeze(cf_pre_P(:,isub,:)),1);
    cf_pre_N_sub_avg = nanmean(squeeze(cf_pre_N(:,isub,:)),1);
    cf_pre_comb_sub_avg = cf_pre_P_sub_avg/2 + -cf_pre_N_sub_avg/2;
    
    ar_P_sub_avg = cf_post_P_sub_avg - cf_pre_P_sub_avg;
    ar_N_sub_avg = cf_post_N_sub_avg - cf_pre_N_sub_avg;
    ar_comb_sub_avg = ar_P_sub_avg/2 + -ar_N_sub_avg/2;
    
    %calculate the standard error for each trial type
    cf_pert_P_sub_err = nanstd(squeeze(cf_pert_P(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_pert_P(:,isub,:))),2))) ;
    cf_pert_N_sub_err = nanstd(squeeze(cf_pert_N(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_pert_N(:,isub,:))),2)));
    cf_pert_comb_sub_err = nanstd([squeeze(cf_pert_P(:,isub,:)); -squeeze(cf_pert_N(:,isub,:))],0,1)/...
        sqrt(sum(~all(isnan([squeeze(cf_pert_P(:,isub,:)); -squeeze(cf_pert_N(:,isub,:))]),2)));
    
    cf_post_P_sub_err = nanstd(squeeze(cf_post_P(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_post_P(:,isub,:))),2)));
    cf_post_N_sub_err = nanstd(squeeze(cf_post_N(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_post_N(:,isub,:))),2)));
    cf_post_comb_sub_err= nanstd([squeeze(cf_post_P(:,isub,:)); -squeeze(cf_post_N(:,isub,:))],0,1)/...
        sqrt(sum(~all(isnan([squeeze(cf_post_P(:,isub,:)); -squeeze(cf_post_N(:,isub,:))]),2)));
    
    cf_pre_P_sub_err = nanstd(squeeze(cf_pre_P(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_pre_P(:,isub,:))),2)));
    cf_pre_N_sub_err = nanstd(squeeze(cf_pre_N(:,isub,:)),0,1)/sqrt(sum(~all(isnan(squeeze(cf_pre_N(:,isub,:))),2)));
    cf_pre_comb_sub_err= nanstd([squeeze(cf_pre_P(:,isub,:)); -squeeze(cf_pre_N(:,isub,:))],0,1)/...
        sqrt(sum(~all(isnan([squeeze(cf_pre_P(:,isub,:)); -squeeze(cf_pre_N(:,isub,:))]),2)));
    
    ar_P_sub_tmp = squeeze(cf_post_P(:,isub,:)) - squeeze(cf_pre_P(:,isub,:));
    ar_N_sub_tmp = squeeze(cf_post_N(:,isub,:)) - squeeze(cf_pre_N(:,isub,:));
    ar_comb_sub_tmp = [ar_P_sub_tmp; -ar_N_sub_tmp];
    ar_comb_sub_err = nanstd(ar_comb_sub_tmp,0,1)/sqrt(sum(~all(isnan(ar_comb_sub_tmp),2)));    
    
    %determine the adaptive response and ideal forces again based on MEDIANS
    cf_pert_P_sub_med = nanmedian(squeeze(cf_pert_P(:,isub,:)),1);
    cf_pert_N_sub_med = nanmedian(squeeze(cf_pert_N(:,isub,:)),1);
    cf_pert_comb_sub_med = cf_pert_P_sub_med/2 + -cf_pert_N_sub_med/2;
    
    cf_post_P_sub_med = nanmedian(squeeze(cf_post_P(:,isub,:)),1);
    cf_post_N_sub_med = nanmedian(squeeze(cf_post_N(:,isub,:)),1);
    cf_post_comb_sub_med = cf_post_P_sub_med/2 + -cf_post_N_sub_med/2;
    
    cf_pre_P_sub_med = nanmedian(squeeze(cf_pre_P(:,isub,:)),1);
    cf_pre_N_sub_med = nanmedian(squeeze(cf_pre_N(:,isub,:)),1);
    cf_pre_comb_sub_med = cf_pre_P_sub_med/2 + -cf_pre_N_sub_med/2;
    
    ar_P_sub_med = cf_post_P_sub_med - cf_pre_P_sub_med;
    ar_N_sub_med = cf_post_N_sub_med - cf_pre_N_sub_med;
    ar_comb_sub_med = ar_P_sub_med/2 + -ar_N_sub_med/2;
    
    %save the data
    ar_sub_all(:,isub) = ar_comb_sub_avg;
    ideal_avg_sub_all(:,isub) = cf_pert_comb_sub_avg;
    ideal_median_sub_all(:,isub) = cf_pert_comb_sub_med;
    pre_sub_all(:,isub) = cf_pre_comb_sub_avg;
    post_sub_all(:,isub) = cf_post_comb_sub_avg;
    
    %determine the learning rate for each subject
    %first do it by regression
    [b_diff_reg_sub, ~] = regress(ar_comb_sub_avg(:), [cf_pert_comb_sub_avg(:), cf_pert_comb_sub_avg(:)*0+1]);
    
    %save them
    LR_reg_sub_all(isub) = b_diff_reg_sub(1);
    
    %plot data
    subplot(2, num_plots, cc); hold on;
    
    %plot the perturbation
    plot(tt, cf_pert_comb_sub_avg, 'color', grey);
    
    %post
    plot(tt, cf_post_comb_sub_avg, 'color', 'r');
    
    %pre
    plot(tt, cf_pre_comb_sub_avg, 'color', 'b');
    
    %plot the adaptive response
    plot(tt, ar_comb_sub_avg, 'g', 'linewidth', 1);
    
    if show_error
        display_xy_error(tt, cf_pert_comb_sub_avg, [], cf_pert_comb_sub_err, 'color', grey);
        display_xy_error(tt, cf_post_comb_sub_avg, [], cf_post_comb_sub_err, 'color', 'r');
        display_xy_error(tt, cf_pre_comb_sub_avg, [], cf_pre_comb_sub_err, 'color', 'b');
        display_xy_error(tt, ar_comb_sub_avg, [], ar_comb_sub_err, 'color', 'g');
    end
    
    %optionally, plot the positive and negative perturbations and flipping
    %the negative one!
    if show_PN
        plot(tt, cf_pert_P_sub_avg, 'color', cyan, 'Linestyle', '-');
        plot(tt, flip_sign*cf_pert_N_sub_avg, 'color', pink, 'Linestyle', '-');
        
        if show_error
            display_xy_error(tt, cf_pert_P_sub_avg, [], cf_pert_P_sub_err, 'color', cyan);
            display_xy_error(tt, flip_sign*cf_pert_N_sub_avg, [], cf_pert_N_sub_err, 'color', pink);
        end
    end
    
    %optionally show medians
    if show_median
        plot(tt, cf_pert_comb_sub_med, 'color', grey*col_gain, 'Linestyle', '-'); %perturbation
        plot(tt, cf_post_comb_sub_med, 'color', [1,0,0]*col_gain, 'Linestyle', '-'); %post
        plot(tt, cf_pre_comb_sub_med, 'color', [0,0,1]*col_gain, 'Linestyle', '-'); %pre
        plot(tt, ar_comb_sub_med, 'color',[0,1,0]*col_gain, 'linewidth', 1, 'Linestyle', '-'); %AR
    end
    
    if show_PN && show_median
        plot(tt, cf_pert_P_sub_med, 'color', cyan*col_gain, 'Linestyle', '-');
        plot(tt, flip_sign*cf_pert_N_sub_med, 'color', pink*col_gain, 'Linestyle', '-');
    end
    
    %plot the grand mean
    %plot(tt, ar_diff_gm_tmp, 'k', 'linewidth', 2.5);
    
    %make the legend
    %         legend('AutoUpdate','Off');
    %         leg_tmp = legend({'perturbation', 'post', 'pre', 'AR'});
    %         set(leg_tmp, 'Location', 'Best');
    
    title([sub_id{isub}, ', LR = ', num2str(LR_reg_sub_all(isub))]);
    %title([sub_id{isub}]);
    
    xlabel('time (sec)');
    ylabel('Force (N)');
    %ylim([force_LB, force_UB]);
    grid on;
    
    cc = cc+1;
end

%%%see if calculating adaptive response from grand mean is much different 
[b_diff_reg_gm, ~] = regress(ar_gm(:), [ideal_gm(:), ideal_gm(:)*0+1]);

%plot the population AR
subplot(2,num_plots, cc); hold on;

plot(nan, nan, 'color', 'g', 'displayname', 'AR');
plot(nan, nan, 'color', grey, 'displayname', 'Ideal');
plot(nan, nan, 'color', 'b', 'displayname', 'pre');
plot(nan, nan, 'color', 'r', 'displayname', 'post');
if show_error
    plot(nan, nan, 'color', cyan, 'displayname', 'positive trials');
    plot(nan, nan, 'color', pink, 'displayname', 'negative trials');
end
% plot(nan, nan, 'color', 'k', 'linestyle', '-', 'displayname', 'Mean');
% plot(nan, nan, 'color', 'k', 'linestyle', '--', 'displayname', 'Median');

plot(tt, ar_gm, 'g', 'linewidth', 1.5);
plot(tt, ideal_gm, 'color', grey, 'linewidth', 1.5);
plot(tt, pre_gm, 'color', 'b', 'linewidth', 1.5);
plot(tt, post_gm, 'color', 'r', 'linewidth', 1.5);

if show_error
    display_xy_error(tt, ar_gm, [], ar_comb_gm_err, 'color', 'g');
    display_xy_error(tt, ideal_gm, [], ideal_gm_err, 'color', grey);
    display_xy_error(tt, pre_gm, [], pre_gm_err, 'color', 'b');
    display_xy_error(tt, post_gm, [], post_gm_err, 'color', 'r');
end

if show_PN
    plot(tt, ideal_P_gm_tmp, 'color', cyan, 'Linestyle', '-');
    plot(tt, flip_sign*ideal_N_gm_tmp, 'color', pink, 'Linestyle', '-');
    
    if show_error
        display_xy_error(tt, ideal_P_gm_tmp, [], ideal_P_gm_err, 'color', cyan);
        display_xy_error(tt, flip_sign*ideal_N_gm_tmp, [], ideal_N_gm_err, 'color', pink);        
    end
end

% if show_median
%     plot(tt, ideal_med_gm, 'color', grey, 'Linestyle', '--'); %perturbation
%     plot(tt, post_med_gm, 'color', 'r', 'Linestyle', '--'); %post
%     plot(tt, pre_med_gm, 'color', 'b', 'Linestyle', '--'); %pre
%     plot(tt, ar_med_gm, 'color','k', 'linewidth', 1, 'Linestyle', '--'); %AR
% end

xlabel('time (sec)');
ylabel('Force (N)');
%ylim([force_LB, force_UB]);
grid on;

legend('AutoUpdate','Off');
if show_error
    leg2_tmp = legend({'AR', 'combined ideal', 'pre', 'post', 'positive ideal trials', 'negative ideal trials'});
else
    leg2_tmp = legend({'AR', 'combined ideal', 'pre', 'post'});
end
set(leg2_tmp, 'Location', 'Best');

%title(['Grand means, LR-reg = ', num2str(nanmean(LR_reg_sub_all)), ' vs ', num2str(b_diff_reg_gm(1))]);
title(['Grand means, LR-reg = ', num2str(nanmean(LR_reg_sub_all)), ' +/- ', num2str(std(LR_reg_sub_all))]);
%title(['Grand means for ',pert_fld]);

h_tmp = suptitle(['Force profiles for ', pert_fld, ' trials in a ', env_fld]);
set(h_tmp, 'Interpreter', 'None');
%keyboard;
end