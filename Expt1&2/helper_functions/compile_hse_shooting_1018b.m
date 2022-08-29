clear all;
close all;
home;

% cd('C:\Users\Laith Alhussein\Desktop\Research\HSE_exp\expb\data\mat_files\ALL_data');
cd('C:\Users\Laith\Dropbox (HNL)\Laith_files\HSE_exp\expb\data\mat_files\ALL_DATA');

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
hse_data=load(strcat(Cdr,'\',hse_filename));
dat=hse_data.dat;
info=hse_data.info;

[total_trials, num_subjects] = size(dat.n);
win = [-79:0];
num_samples = length(win);
original_samples = info.num_samples;
B = info.FFMAG; % 7.5 as of 8/27/2018
tt = [0:length(win)-1]*0.005;
EC_err_size = 9;

%colors
purple = [0.5,0,0.5];
grey = [0.5,0.5,0.5];

%% get the correct indices for each participant

%specifiy the blocks for each period, including both directions for now
% fam_block = 1;
% baseline_block = [2,3];
% training1_block = [4,5];
% test1_block = [6,7,8,9,10,11];
% test2_block = [12,13,14,15,16,17];
% training3_block = [18,19];

fam_idx = [1:60];
baseline_idx = [fam_idx(end)+ 1 : fam_idx(end) + 200];
training1_idx = [baseline_idx(end) + 1 : baseline_idx(end) + 300];
test1_idx = [training1_idx(end) + 1 : training1_idx(end) +  990];
test2_idx = [test1_idx(end) + 1 : test1_idx(end) + 990];
training2_idx = [test2_idx(end) + 1: test2_idx(end) + 300];

%concentrate on forward trials...
fam_idx = fam_idx(2:2:end)/2;
baseline_idx = baseline_idx(2:2:end)/2;
training1_idx = training1_idx(2:2:end)/2;
test1_idx = test1_idx(2:2:end)/2;
test2_idx = test2_idx(2:2:end)/2;
training2_idx = training2_idx(2:2:end)/2;

tpb = [length(fam_idx), length(baseline_idx), length(training1_idx), length(test1_idx), length(test2_idx),...
    length(training2_idx)];
exp_seq_all = cumsum(tpb);
exp_seq = cumsum(tpb(3:end));

tgt_all = cell(num_subjects,1);
fwrd_tgt_all = cell(num_subjects,1);

for k1=1:num_subjects
    
    idx.baseline.all(k1,:) = baseline_idx;
    
    idx.training1(k1,:) = training1_idx;
    idx.test1.all(k1,:) = test1_idx;
    
    idx.training2(k1,:) = training2_idx;
    idx.test2.all(k1,:) = test2_idx;
    
    sub_tgt_tmp = info.tgt_dat(k1,:);
    sub_tgt_all = vertcat(sub_tgt_tmp{:});
    sub_tgt_all_fwrd = sub_tgt_all(2:2:end,:);
    
    fwrd_tgt_all{k1} = sub_tgt_all_fwrd;
    tgt_all{k1} = sub_tgt_all_fwrd;
    
    ec_tmp1 = find(sub_tgt_all_fwrd(:,7) == 1);
    %ec_tmp2 = find(sub_tgt_all_fwrd(:,7) == 2 & abs(sub_tgt_all_fwrd(:,end))~=EC_err_size ); %not counting the vEC triplet
    ec_tmp2 = find(sub_tgt_all_fwrd(:,7) == 2 ); %not counting the vEC triplet
    
    idx.baseline.ec(k1,:) = ec_tmp1(ec_tmp1 <= idx.baseline.all(k1,end));
    
    idx.vEC.test1(k1,:) = ec_tmp2( ec_tmp2<=idx.test1.all(k1,end) & ec_tmp2>=idx.test1.all(k1,1) );
    idx.vEC.test2(k1,:) = ec_tmp2( ec_tmp2<=idx.test2.all(k1,end) & ec_tmp2>=idx.test2.all(k1,1) );
    
    %find the triplets for each case
    triplets_tmp1 = find(sub_tgt_all_fwrd(:,2)~=0);
    triplets_tmp2 = find(abs(sub_tgt_all_fwrd(:,end))==EC_err_size);
    
    idx.err_FF.all = triplets_tmp1;
    idx.err_EC.all = triplets_tmp2;
    
    idx.pre.err_FF(k1,:) = triplets_tmp1 -1;     idx.post.err_FF(k1,:) = triplets_tmp1+1;
    idx.pre.err_EC(k1,:) = triplets_tmp2 -1;     idx.post.err_EC(k1,:) = triplets_tmp2 +1;
    
    %find the sign of the applied errors during the triplets
    err_sign.EC(k1,:) = sub_tgt_all_fwrd(triplets_tmp1,2);
    err_sign.FF(k1,:) = sign(sub_tgt_all_fwrd(triplets_tmp2,end));
    
    %combine all vEC trials together
    idx.vEC.all(k1,:) = ec_tmp2;
    
    %get all non-baseline EC trials
    idx.EC.all(k1,:) = ec_tmp1(ec_tmp1 > idx.baseline.all(k1,end));
    
end

%% align data to the target

vel = cell(total_trials/2, 1);
ypos = cell(total_trials/2, 1);
xpos = cell(total_trials/2, 1);
fr = cell(total_trials/2, 1);

err = nan(num_subjects, total_trials/2); %we should also calculate the error that the participant experienced
%vmax = err;

short_movements = 0;
junk_forces = 0;
for k2=1:num_subjects
    cc=1;
    for k3=2:2:total_trials
        
        tgt_idx = find( dat.pxr{k3}(:,k2) > 0.1, 1, 'first') - 1;
        
        force_tmp = dat.fyr{k3};
        
        %we can have cases in the early trials where there is no force data saved (these are just null trials anyway)
        if size(force_tmp,2)~=num_subjects, force_tmp(:,end:num_subjects) = NaN; end
        
        %sometimes the movement can be short, and we dont have enough data
        %before the pillow to subsample the movement
        if ~isempty(tgt_idx) & tgt_idx<length(win)
            vel{cc}(:,k2) = nan+win;
            ypos{cc}(:,k2) = nan+win;
            xpos{cc}(:,k2) = nan+win;
            if ~isempty(force_tmp)
                fr{cc}(:,k2) = nan+win;
            else
                fr{cc} = [];
            end
            short_movements = short_movements+1;
            
        else
            if isempty(tgt_idx), tgt_idx = abs(min(win)) + 1; end
            
            vel{cc}(:,k2) = dat.vxr{k3}(tgt_idx+win,k2);
            ypos{cc}(:,k2) = dat.pxr{k3}(tgt_idx+win,k2);
            xpos{cc}(:,k2) = dat.pyr{k3}(tgt_idx+win,k2);
            
            err(k2, cc) = dat.pyr{k3}(tgt_idx,k2) * 1000; %keep in mm
            vmax.all(k2,cc) = dat.vmax(k3,k2);

            if isempty(force_tmp)
                fr{cc}= [];
            else
                fr{cc}(:,k2) = force_tmp(win+tgt_idx, k2);
            end
            
        end
        
        cc=cc+1;
    end
end

%% now calculate the adaptation on all trials by integration

ILF_all = nan(num_subjects, total_trials/2);

for k6=1:total_trials/2
    if ~isempty(fr{k6})
        for k7=1:num_subjects
            %ACI_all(k7,k6) = trapz(fr{k6}(:,k7)) / trapz(vel{k6}(:,k7) * B);
            ILF_all(k7,k6) = sum(fr{k6}(:,k7)) / sum( vel{k6}(:,k7) * B);
        end
    end
end

%% subtract baseline AC's

baseline_subtract_flag = 1;

for kq = 1:num_subjects
    
    ILF.baseline(kq,:) = ILF_all(kq, idx.baseline.ec(kq,:));
    
    if baseline_subtract_flag 
        %keyboard;
        ILF_all(kq,:) = ILF_all(kq,:) - nanmean(ILF.baseline(kq,:));        
        baseline_subtract_flag = 0;
    end
    
end

%% subtract each participant's mean response over all EC (vEC and zEC) trials

mean_subtract_flag = 1;

if mean_subtract_flag
    
    %ILF_all = bsxfun(@(x,y) x-y, ILF_training1_cut_raw, nanmean(ILF_training1_cut_raw,2));
    
    for k=1:num_subjects
        
        ILF_all(k, [idx.vEC.all(k,:), idx.EC.all(k,:)]) = ILF_all(k, [idx.vEC.all(k,:), idx.EC.all(k,:)]) - ...
            nanmean(ILF_all(k, [idx.vEC.all(k,:), idx.EC.all(k,:)]) );
    end
    mean_subtract_flag = 0;
end

%% save the AC data in a struct

ILF_unflipped_all = ILF_all;

for kq=1:num_subjects
    
    ILF.pre.err_FF(kq,:) = ILF_all(kq, idx.pre.err_FF(kq,:));
    ILF.pre.err_EC(kq,:) = ILF_all(kq, idx.pre.err_EC(kq,:));
    
    ILF.post.err_FF(kq,:) = ILF_all(kq, idx.post.err_FF(kq,:));
    ILF.post.err_EC(kq,:) = ILF_all(kq, idx.post.err_EC(kq,:));
    
    ILF.ar.err_FF(kq,:) = ( ILF.post.err_FF(kq,:) - ILF.pre.err_FF(kq,:) ) .* -err_sign.FF(kq,:);
    ILF.ar.err_EC(kq,:) = ( ILF.post.err_EC(kq,:) - ILF.pre.err_EC(kq,:) ) .* -err_sign.EC(kq,:);
    
    %save the adaptive response without flipping it
    ILF_uf.pre.err_FF(kq,:) = ILF_all(kq, idx.pre.err_FF(kq,:));
    ILF_uf.pre.err_EC(kq,:) = ILF_all(kq, idx.pre.err_EC(kq,:));
    
    ILF_uf.post.err_FF(kq,:) = ILF_all(kq, idx.post.err_FF(kq,:));
    ILF_uf.post.err_EC(kq,:) = ILF_all(kq, idx.post.err_EC(kq,:));
    
    ILF_uf.ar.err_FF(kq,:) = ( ILF_uf.post.err_FF(kq,:) - ILF_uf.pre.err_FF(kq,:) );
    ILF_uf.ar.err_EC(kq,:) = ( ILF_uf.post.err_EC(kq,:) - ILF_uf.pre.err_EC(kq,:) );
    
    %save the training data
    ILF.training1(kq,:) = ILF_all(kq,idx.training1(kq,:));
    ILF.training2(kq,:) = ILF_all(kq,idx.training2(kq,:));
    
    %test data as well, with ALL trials
    ILF.test1(kq,:) = ILF_all(kq,idx.test1.all(kq,:));
    ILF.test2(kq,:) = ILF_all(kq,idx.test2.all(kq,:));
    
end

%% make a time series plot of the vEC trials only, but make the discontinuities apparent

figure; hold on;

%get all trials from training onwards
%ILF_all_tmp = ILF_all(:, training1_idx(1):end);
ILF_all_tmp = ILF_all;

%replace the triplet trials with nan
ILF_vEC_all = ILF_all_tmp;
ILF_vEC_all(:, [idx.err_FF.all, idx.err_FF.all+1, idx.err_FF.all-1,...
    idx.err_EC.all, idx.err_EC.all+1, idx.err_EC.all-1]) = NaN;

%take it only from training onwards
ILF_vEC_avg = nanmean(ILF_vEC_all(:, training1_idx(1):end), 1);

plot(ILF_vEC_avg);

ylabel('ILF');
xlabel('Trial # wrt training onset');

yy1 = [-6,5];

ylim(yy1);

plot(ones(1,2).*exp_seq(1),yy1, 'k--'); %end of first training period
plot(ones(1,2).*exp_seq(2),yy1, 'k--');
plot(ones(1,2).*exp_seq(3),yy1, 'k--');


%% plot the pre and post trials
figure;
LW=2;
num_trip = 33;
x1 = [1:num_trip];

for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(ILF.pre.err_FF(k,:), 'k', 'displayname', 'FF');  plot(x1,mean(ILF.pre.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF.pre.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF.pre.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(ILF_uf.pre.err_FF(k,:), 'k', 'displayname', 'FF'); plot(x1,mean(ILF_uf.pre.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF_uf.pre.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF_uf.pre.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);

end

suptitle('Pre AC trials');

figure;
for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(ILF.post.err_FF(k,:), 'k', 'displayname', 'FF');  plot(x1,mean(ILF.post.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF.post.err_EC(k,:), 'color',grey, 'displayname', 'EC'); plot(x1,mean(ILF.post.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(ILF_uf.post.err_FF(k,:), 'k', 'displayname', 'FF'); plot(x1,mean(ILF_uf.post.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF_uf.post.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF_uf.post.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
end

suptitle('Post AC trials');

%%
figure;
for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(ILF.ar.err_FF(k,:), 'k', 'displayname', 'FF');  plot(x1,mean(ILF.ar.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF.ar.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF.ar.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(ILF_uf.ar.err_FF(k,:), 'k', 'displayname', 'FF'); plot(x1,mean(ILF_uf.ar.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF_uf.ar.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF_uf.ar.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);

    ylim([-4 4]);
    
end

for k=1:num_subjects,
    subplot(2,num_subjects,k); hold on;
    plot(ILF.ar.err_FF(k,:), 'k', 'displayname', 'FF');  plot(x1,mean(ILF.ar.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF.ar.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF.ar.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
    
    subplot(2,num_subjects,k+num_subjects); hold on;
    plot(ILF_uf.ar.err_FF(k,:), 'k', 'displayname', 'FF'); plot(x1,mean(ILF_uf.ar.err_FF(k,:))*ones(1,num_trip), 'k','linewidth',2);
    plot(ILF_uf.ar.err_EC(k,:), 'color', grey, 'displayname', 'EC'); plot(x1,mean(ILF_uf.ar.err_EC(k,:))*ones(1,num_trip), 'color',grey,'linewidth',2);
    
    xlabel([info.sublist{k}]);
end

for k=1:num_subjects*2, subplot(2,num_subjects,k); ylim([-3 3]); end

suptitle('Adaptive response');

%% plot the force profiles, with the ideal, for each participant

cf1 = 'err_FF';

for i=1:num_subjects,
    figure; 
    for k=1:num_trip
        
        subplot(ceil(num_trip/num_subjects),num_subjects,k);hold on;
        
        %plot the actual pre and post forces
        
        plot(fr{idx.pre.(cf1)(i,k)}(:,i), 'r');
        plot(fr{idx.post.(cf1)(i,k)}(:,i), 'b');
        
        %plot the ideal
        plot(vel{idx.pre.(cf1)(i,k)}(:,i) * B *-err_sign.(cf1(end-1:end))(i,k) , 'k');
        
        %plot the AR
        plot(fr{idx.post.(cf1)(i,k)}(:,i) - fr{idx.pre.(cf1)(i,k)}(:,i), 'color', [0.5, 0, 0.5]);
        
        grid on;
        
        title([cf1, ', ILF is: ', num2str(ILF.ar.(cf1)(i,k),2)] );
        
    end
    suptitle(['Pre & Post: ', info.sublist{i}]);
end                            

%% calculate the average FP for each case

f_sub.pre.err_FF = get_sub_data(fr, idx.pre.err_FF);
f_sub.pre.err_EC = get_sub_data(fr, idx.pre.err_EC);

f_sub.post.err_FF = get_sub_data(fr, idx.post.err_FF);
f_sub.post.err_EC = get_sub_data(fr, idx.post.err_EC);

f_sub.ar.err_FF = f_sub.post.err_FF - f_sub.pre.err_FF;
f_sub.ar.err_EC = f_sub.post.err_EC - f_sub.pre.err_EC;

vel_sub.pre.err_FF = get_sub_data(vel, idx.pre.err_FF);
vel_sub.pre.err_EC = get_sub_data(vel, idx.pre.err_EC);

vel_sub.post.err_FF = get_sub_data(vel, idx.post.err_FF);
vel_sub.post.err_EC = get_sub_data(vel, idx.post.err_EC);

ypos_sub.pre.err_FF = get_sub_data(ypos, idx.pre.err_FF);
ypos_sub.pre.err_EC = get_sub_data(ypos, idx.pre.err_EC);

ypos_sub.post.err_FF = get_sub_data(ypos, idx.post.err_FF);
ypos_sub.post.err_EC = get_sub_data(ypos, idx.post.err_EC);


%% flip the sign of the FF for the adaptive response

f2_sub.ar.err_FF = f_sub.ar.err_FF * NaN;
f2_sub.ar.err_EC = f_sub.ar.err_EC * NaN;

for q1=1:num_subjects
    
    %if q1==4, keyboard; end
    
    for q2 = 1:num_trip
        f2_sub.ar.err_FF(q1,q2,:) = f_sub.ar.err_FF(q1,q2,:) * -err_sign.FF(q1,q2);
        f2_sub.ar.err_EC(q1,q2,:) = f_sub.ar.err_EC(q1,q2,:) * -err_sign.EC(q1,q2); 
    end
    
end

%% instead of including pre and post, just make 1 plot with all subjects and plot the ar
%make sure to add the AC...

num_early_trials = 1;
flds = {'err_FF', 'err_EC'};
num_types = length(flds);

figure;

for p1 = 1:num_types
    for kk=1:num_subjects
        
        cfld = flds{p1};
        
        %first do vEC
        subplot(num_types,num_subjects, kk + (p1-1) * num_subjects); hold on;
        for jj=1:num_early_trials
            plot(squeeze(smooth( f2_sub.ar.(cfld)(kk,jj,:), 7) ) , 'color', purple);
        end
        
        %show the mean AR and ideal (get ideal from pre trials)
        plot( squeeze( nanmean(f2_sub.ar.(cfld)(kk, [1:num_early_trials], :), 2) ), 'color', purple, 'linewidth', 2);
        plot( squeeze( nanmean(vel_sub.pre.(cfld)(kk, [1:num_early_trials], :) * B, 2) ), 'k', 'linewidth', 2);
        
        ylim([-10 15]); grid on;
        
        if p1==1, title([info.sublist{kk}]); end
        
        if kk==1
            if p1==1
                ylabel(cfld(end-1:end));
            elseif p1==2
                ylabel(cfld(end-1:end));
            end
        end
    end
end

suptitle('Early Adaptive Response');

%keyboard;


%% first look to see how well the imposed errors correlate with the experienced errors

%get the imposed error sequence in mm (separate based on training only or training + test)
ierr_tmp = tgt_all{1}(idx.vEC.all(1,:), end); %this is in degrees
ierr.all  = tand(ierr_tmp) *100; %in mm

%repeat ierr since we will combine participant data
ierr_rep = repmat(ierr.all, [num_subjects,1]);

%concatinate the participant data
err_vEC_all = nan(num_subjects * length(ierr_tmp), 1);

c1 = 1;
for k=1:num_subjects
    for j=1:length(ierr_tmp)
        err_vEC_all(c1) = err(k, idx.vEC.all(k,j));
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


%% look at the AC, and do the modeling for both vEC training periods

burn_trials = 40; %a certain number of trials need to be removed to exclude the effects of the EC bias

ILF_training1_cut = ILF.training1(:, [burn_trials+1:end]);
ILF_training2_cut = ILF.training2(:, [burn_trials+1:end]);

ILF_test1 = ILF.test1;
ILF_test2 = ILF.test2;

%get the imposed errors for each case in mm
ierr_training1d_tmp = -tgt_all{1}(idx.training1(1,:), end); %multiply by -1 to let CCW be leftward
ierr_training1_tmp = tand(ierr_training1d_tmp) * 100;
ierr_training1_cut = ierr_training1_tmp([burn_trials+1:end]);

ierr_training2d_tmp = -tgt_all{1}(idx.training2(1,:), end);
ierr_training2_tmp = tand(ierr_training2d_tmp) * 100;
ierr_training2_cut = ierr_training2_tmp([burn_trials+1:end]);

ierr_test1d_tmp = -tgt_all{1}(idx.test1.all(1,:), end);
ierr_test1 = tand(ierr_test1d_tmp) * 100;

ierr_test2d_tmp = -tgt_all{1}(idx.test2.all(1,:), end);
ierr_test2 = tand(ierr_test2d_tmp) * 100;

%calculate subject mean
ILF_avg_training1 = nanmean(ILF_training1_cut, 1);
ILF_avg_training2 = nanmean(ILF_training2_cut, 1);

ILF_avg_test1 = nanmean(ILF_test1,1);
ILF_avg_test2 = nanmean(ILF_test2,1);

%%%%% to estimate the effect of stiffness, regress the adpatation onto the imposed error encountered on the same trial
%training 1
[b_stiff_tr1,b_stiff_int_tr1,res_tr1,~,s_stiff_tr1] = regress(ILF_avg_training1', ierr_training1_cut);

%training 2
[b_stiff_tr2,b_stiff_int_tr2,res_tr2,~,s_stiff_tr2] = regress(ILF_avg_training2', ierr_training2_cut);

%test 1
[b_stiff_tst1, b_stiff_int_tst1, res_tst1,~,s_stiff_tst1] = regress(ILF_avg_test1', ierr_test1);

%test 2
[b_stiff_tst2, b_stiff_int_tst2, res_tst2,~,s_stiff_tst2] = regress(ILF_avg_test2', ierr_test2);

%use the estimate of stiffness over each period to convert the error to normalized errors
ierr_tr1_norm = b_stiff_tr1*ierr_training1_cut;
ierr_tst1_norm = b_stiff_tst1*ierr_test1;
ierr_tst2_norm = b_stiff_tst2*ierr_test2;
ierr_tr2_norm = b_stiff_tr2*ierr_training2_cut;

%plot participant mean
figure; suptitle('ILF (averaged)');

subplot(411); hold on; ylabel('Training 1');
plot(ierr_training1_cut, ILF_avg_training1, '.');

subplot(412); hold on; ylabel('Test 1');
plot(ierr_test1, ILF_avg_test1, '.');

subplot(413); hold on; ylabel('Test 2');
plot(ierr_test2, ILF_avg_test2, '.');

subplot(414); hold on; ylabel('Training 2');
plot(ierr_training2_cut, ILF_avg_training2, '.');

xlabel('Error (mm)');

%axis([-10, 10, -1, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corr(ierr_training1_rep, ILF_training1_cat, 'rows', 'complete')
% corr(ierr_training1_cut, ILF_avg_training1')
% 
% corr(ierr_training2_rep, ILF_training2_cat, 'rows', 'complete')
% corr(ierr_training2_cut, ILF_avg_training2')

%% estimate the learning after removing stiffness

%qq = nlinfit(ierr_training1_cut, ILF_avg_training1, @stiffness_model, [0.8,0.07,0.07]);

%%% fit the model

%training 1
[p_tr1,rr_tr1,J_tr1] = nlinfit(ierr_tr1_norm, res_tr1, @stiffness_model2, [0.8,0.07,0.1]);
ci_tr1 = nlparci(p_tr1,rr_tr1,'Jacobian',J_tr1);

%training 2
[p_tr2,rr_tr2,J_tr2] = nlinfit(ierr_tr2_norm, res_tr2, @stiffness_model2, [0.8,0.07,0.1]);
ci_tr2 = nlparci(p_tr2,rr_tr2,'Jacobian',J_tr2);

%test 1
[p_tst1,rr_tst1,J_tst1] = nlinfit(ierr_tst1_norm, res_tst1, @stiffness_model2, [0.8,0.07,0.1]);
ci_tst1 = nlparci(p_tst1,rr_tst1,'Jacobian',J_tst1);

%calculate the r-squared
tst1_fit = stiffness_model2(p_tst1, ierr_tst1_norm);
corr_tst1 = corrcoef(tst1_fit, res_tst1); 
rsq_tst1 = corr_tst1(1,2)^2;

%test 2
[p_tst2,rr_tst2,J_tst2] = nlinfit(ierr_tst2_norm, res_tst2, @stiffness_model2, [0.8,0.07,0.1]);
ci_tst2 = nlparci(p_tst2,rr_tst2,'Jacobian',J_tst2);

elr = [p_tr1(2), p_tst1(2), p_tst2(2), p_tr2(2)]; 

%fit it onto a model that saturates with error size
%qq3 = nlinfit(ierr_training_cut, ILF_avg_t1_ns', @stiffness_model_sat, [0.8, 0.07, 5, 0.1] );

%% estimate the LR on individual participant basis

%b_stiff.tr1 = [];
LR_all = nan(num_subjects, 4); %4 model paramters

for kq=1:num_subjects
    
    % estimate stiffness first
    [b_stiff.tr1(kq), ~, res_stiff.tr1(kq,:)] = regress( ILF_training1_cut(kq,:)', ierr_training1_cut );
    [b_stiff.tr2(kq), ~, res_stiff.tr2(kq,:)] = regress( ILF_training2_cut(kq,:)', ierr_training2_cut );
    
    [b_stiff.tst1(kq), ~, res_stiff.tst1(kq,:)] = regress( ILF_test1(kq,:)', ierr_test1 );
    [b_stiff.tst2(kq), ~, res_stiff.tst2(kq,:)] = regress( ILF_test2(kq,:)', ierr_test2 );
    
    %calculate the normalized error
    norm_err.tr1(kq,:) = b_stiff.tr1(kq)*ierr_training1_cut;
    norm_err.tr2(kq,:) = b_stiff.tr2(kq)*ierr_training2_cut;
    norm_err.tst1(kq,:) = b_stiff.tst1(kq)*ierr_test1;
    norm_err.tst2(kq,:) = b_stiff.tst2(kq)*ierr_test2;
    
    %estimate the learning rate
    [p_lr.tr1(kq,:)] = nlinfit(norm_err.tr1(kq,:), res_stiff.tr1(kq,:)', @stiffness_model2, [0.8,0.07,0.1]);
    [p_lr.tr2(kq,:)] = nlinfit(norm_err.tr2(kq,:), res_stiff.tr2(kq,:)', @stiffness_model2, [0.8,0.07,0.1]);
    
    [p_lr.tst1(kq,:)] = nlinfit(norm_err.tst1(kq,:), res_stiff.tst1(kq,:)', @stiffness_model2, [0.8,0.07,0.1]);
    [p_lr.tst2(kq,:)] = nlinfit(norm_err.tst2(kq,:), res_stiff.tst2(kq,:)', @stiffness_model2, [0.8,0.07,0.1]);
    
    LR_all(kq,:) = [p_lr.tr1(kq,2), p_lr.tst1(kq,2), p_lr.tst2(kq,2), p_lr.tr2(kq,2)];

end

%keyboard;
%% look at displacement on FF trials

%displacement from FF triplets
err_FF_all_sub = abs(err(:,idx.err_FF.all)); %dont care about direction
err_FF_all = reshape(err_FF_all_sub,[size(err_FF_all_sub,1)*size(err_FF_all_sub,2),1]);

err_FF_sub_avg = nanmean(err_FF_all_sub,2);

%displacement from EC triplets
err_EC_all_tmp = abs(err(:,idx.err_EC.all));
err_EC_all = reshape(err_EC_all_tmp, [size(err_EC_all_tmp,1)*size(err_EC_all_tmp,2),1]);

%displacement from vEC's
%use err_vEC_all

%overlay the histograms
figure; hold on;
create_hist_1012_2018(err_FF_all,'k',20);
create_hist_1012_2018(err_EC_all,grey,15);

title('Distribution of displacements');

%% compare max velocities

% vmax.err_FF = reshape(vmax.all(:,idx.err_FF.all), [num_subjects*num_trip,1]);
% vmax.err_EC = reshape(vmax.all(:,idx.err_EC.all),[num_subjects*num_trip,1]);
% 
% %overlay the histograms
% figure; hold on;
% create_hist_1012_2018(vmax.err_FF,'k',20);
% create_hist_1012_2018(vmax.err_EC,grey,15);
% 
% 
% title('Distribution of max velocity');
% 
% 
% %% plot x-position profile
% 
% xpos_sub.err_FF = get_sub_data(xpos, repmat(idx.err_FF.all',num_subjects,1))*1000;
% xpos_sub.err_EC = get_sub_data(xpos, repmat(idx.err_EC.all',num_subjects,1))*1000;
% 
% xp_FF_pa = nanmean(squeeze(nanmedian(xpos_sub.err_FF,2)),1); %population average
% xp_EC_pa = nanmean(squeeze(nanmedian(xpos_sub.err_EC,2)),1);
% 
% xp_FF_std = nanstd(squeeze(nanmedian(xpos_sub.err_FF,2)),0,1);
% xp_EC_std = nanstd(squeeze(nanmedian(xpos_sub.err_EC,2)),0,1);
% 
% figure; hold on;
% 
% plot(tt,xp_FF_pa,'k','linewidth',2.5);
% plot(tt,xp_EC_pa,'color',grey,'linewidth',2.5);
% 
% standard_error_shading_07_16_2015(xp_FF_pa, xp_FF_std, tt, num_subjects, 'k');
% standard_error_shading_07_16_2015(xp_EC_pa, xp_EC_std, tt, num_subjects, grey);
% 
% ylabel('x-position (mm)');
% xlabel('Time (ms)');

%% on the training and test 1 data, manually check learning rate with the following

%b = (xn - xp)/(ec / ec*) + k(en - ep) + (1-A^2)xp

%try it on an individual participant basis
%ierr_training1_cut
%training 1 data

%ILF_avg_training1

% nn1 = size(res_stiff.tr1,2);
% num_rem1 = rem(nn1,3) + 1;
% 
% b_all_tr1 = nan(num_subjects, nn1-num_rem1);
% b_tr1_cut = b_all_tr1(:,2:end); %106 trials
% 
% iqr_thresh = 0;
% 
% err_cutoff = 0.42;
% 
% norm_err.tr1 = (repmat(ierr_training1_cut,1,num_subjects) / diag(err_FF_sub_avg))'; %normalized error
% 
% figure;
% for k=1:num_subjects
%     for i=2: nn1 - num_rem1
%         %keyboard;
%         ar_tmp = (res_stiff.tr1(k,i+1) - res_stiff.tr1(k,i-1)) / abs(norm_err.tr1(k,i)); %should not be stiffness removed
%         %ar_tmp = (res_stiff.tr1(k,i+1) - res_stiff.tr1(k,i-1)) / (abs(ierr_training1_cut(i))/err_FF_sub_avg(k));
%         stiffness_effect_tmp = b_stiff.tr1(k) * (ierr_training1_cut(i+1) - ierr_training1_cut(i-1) );
%         forgetting_effect_tmp = (1 - (p_lr.tr1(k,1)^2))*res_stiff.tr1(k,i-1);
%         
%         lr_tmp = ar_tmp+stiffness_effect_tmp+forgetting_effect_tmp;
%         
%         b_all_tr1(k,i) = lr_tmp;
%         
%         %b_all_tr1(k,i)
%         
%         %if k==1 & i==17, keyboard; end
%                
%     end
%     
%     %estimate the learning rate in 2 ways: (1) median and (2) remove data x iqr from the mean, and then take the mle
%     
%     b_tmp = b_all_tr1(k,:);
%     b_tr1_cut(k,:) = b_all_tr1(k,2:end);
%     
%     %calculate cutoff based on iqr
%     if iqr_thresh>0
%         b_tmp( (b_tmp< nanmean(b_tmp) - iqr_thresh*iqr(b_tmp)) | ((b_tmp> nanmean(b_tmp) + iqr_thresh*iqr(b_tmp))) ) = NaN;
%     else
%         b_tmp = b_all_tr1(k,:);
%     end
%     mle_estimate1 = mle(b_tmp);
%     mu1 = mle_estimate1(1);
%     
%     subplot(1,num_subjects,k); hold on;
%     
%     plot(norm_err.tr1(k,2:end-num_rem1),b_all_tr1(k,[2:end]),'.'); %plot by normalized error
%     %plot(ones(size(b_all_tr1,2),1)*nanmedian(b_all_tr1(k,:)),'color', 'k');
%     %plot(ones(size(b_all_tr1,2),1)*mu1,'color', grey);
%     
%     %cutoff the data based on error size
%     
%     cut_idx = find(norm_err.tr1(k,2:end-num_rem1)>-err_cutoff & norm_err.tr1(k,2:end-num_rem1)<err_cutoff);
%     b_tr1_cut(k,cut_idx) = NaN;
%     
%     %keyboard;
%     %     legend( {'Individual Trials', ['median = ', num2str(nanmedian(b_all_tr1(k,:)))], ['MLE is ', num2str(mu1)] });
%         title(['Learning rate from model is: ', num2str(LR_all(k,1)) ]);
%     
%     plot(norm_err.tr1(k,2:end-num_rem1), b_tr1_cut(k,:), 'o');
%     
%     ylabel('Learning rate');
%     xlabel('normalized error (mm)');
%     
% end
% suptitle('Training 1 data');


%% repeat the above, but for the test1 dataset

nn2 = size(res_stiff.tst1,2);
num_rem2 = rem(nn2,3) + 1;

b_all_tst1 = nan(num_subjects, nn2 - num_rem2);
b_tst1_cut = b_all_tst1(:,2:end);

%norm_err.tst1 = (repmat(ierr_test1,1,num_subjects) / diag(err_FF_sub_avg))'; %normalized error

err_cutoff = 0.42;

EC_err_mm = tand(EC_err_size)*100;

figure;
for k=1:num_subjects
    for i=2: nn2 - num_rem2
        
        %ar_tmp = (res_stiff.tst1(k,i+1) - res_stiff.tst1(k,i-1)) / abs(norm_err.tst1(k,i)); %ILF_avg_training1
        ar_tmp = (ILF_test1(k,i+1) - ILF_test1(k,i-1)) / abs(ierr_tst1_norm(i));
        stiffness_effect_tmp = b_stiff.tst1(k) * (ierr_test1(i+1) - ierr_test1(i-1) );
        forgetting_effect_tmp = (1 - (p_lr.tst1(k,1)^2))*ILF_test1(k,i-1);
        
        lr_tmp = ar_tmp+stiffness_effect_tmp+forgetting_effect_tmp;
        
        b_all_tst1(k,i) = lr_tmp;
        
        %b_all_tr1(k,i)
        
        %if k==1 & i==17, keyboard; end
        
    end
    
    %estimate the learning rate in 2 ways: (1) median and (2) remove data x iqr from the mean, and then take the mle
    
    b_tmp = b_all_tst1(k,:);
    b_tst1_cut(k,:) = b_all_tst1(k,2:end);
    
    subplot(1,num_subjects,k); hold on;
    
    plot(ierr_tst1_norm(2:end-num_rem2),b_all_tst1(k,[2:end]),'.'); %plot by normalized error
    %plot(ones(size(b_all_tst1,2),1)*nanmedian(b_all_tst1(k,:)),'color', 'k');
    %plot(ones(size(b_all_tst1,2),1)*mu1,'color', grey);
    
    %legend( {'Individual Trials', ['median = ', num2str(nanmedian(b_all_tst1(k,:)))], ['MLE is ', num2str(mu1)] });
    title(['Learning rate from model is: ', num2str(LR_all(k,2)) ]);
    
    cut_idx = find(ierr_tst1_norm(2:end-num_rem2)>-err_cutoff & ierr_tst1_norm(2:end-num_rem2)<err_cutoff);
    b_tst1_cut(k,cut_idx) = NaN;
%     
    plot(ierr_tst1_norm(2:end-num_rem2), b_tst1_cut(k,:), 'o');
    
    %add a dashed line for the triplet error
    %plot([EC_err_mm/err_FF_sub_avg(k), EC_err_mm/err_FF_sub_avg(k)], [min(b_tmp(b_tmp~=-Inf)), max(b_tmp(b_tmp~=Inf))], 'k--', 'linewidth', 0.5);
    %plot(-1*[EC_err_mm/err_FF_sub_avg(k), EC_err_mm/err_FF_sub_avg(k)], [min(b_tmp(b_tmp~=-Inf)), max(b_tmp(b_tmp~=Inf))], 'k--', 'linewidth', 0.5);
    plot([EC_err_mm*b_stiff_tst1, EC_err_mm*b_stiff_tst1], [min(b_tst1_cut(k,:)), max(b_tst1_cut(k,:))], 'k--', 'linewidth', 0.5);
    plot(-1*[EC_err_mm*b_stiff_tst1, EC_err_mm*b_stiff_tst1], [min(b_tst1_cut(k,:)), max(b_tst1_cut(k,:))], 'k--', 'linewidth', 0.5);
    
    ylabel('Learning rate');
    xlabel('normalized error (mm)');
    
end

suptitle('Test 1 data');

%%%% calculate average learning rate for trials with the EC error

EC_err_idx = find( abs(ierr_test1(2:end-num_rem2)) == EC_err_mm);
%keyboard;
%nanmedian(b_tst1_cut(:,EC_err_idx),2)

%% calculate the decay from the pre trials

for i=1:num_subjects
    
   ILF.decay.err_EC(i,:) = (1 - (p_lr.tst1(k,1)^2)) * ILF.pre.err_EC(i,:);
   %ILF_uf.decay.err_EC(i,:) = (1 - (p_lr.tst1(k,1)^2)) * ILF_uf.pre.err_EC(i,:);

end

%%

%plot the learning rate from here on the same plot as the one with the post-pre

decay_flag = 1;

figure;
for k=1:num_subjects
    
    subplot(1,num_subjects, k); hold on;
    
    %%%% results from estimating the AR
    plot(x1, b_tst1_cut(k,EC_err_idx),'k', 'linewidth', 1.5);
    
    %show the mean
    plot(x1*0+1*nanmean(b_tst1_cut(k,EC_err_idx) ),'k');
    
    %%%%% results from experimentally measuring the AR
    
    %plot(ILF.ar.err_EC(k,:), 'color', grey);
    plot(x1, ILF_uf.ar.err_EC(k,:), 'color', grey);
    
    %plot(x1*0+1*nanmean(ILF.ar.err_EC(k,:)), 'color', grey);
    plot(x1*0+1*nanmean(ILF_uf.ar.err_EC(k,:)), 'color', grey);
    
    if decay_flag
        %plot(ILF.ar.err_EC(k,:) + ILF.decay.err_EC(k,:), 'color', grey,'linestyle', '--');
        plot(x1, ILF_uf.ar.err_EC(k,:) + ILF.decay.err_EC(k,:), 'color', grey,'linestyle', '--');
        
        %plot(x1*0+1*nanmean(ILF.ar.err_EC(k,:)+ ILF.decay.err_EC(k,:)), 'color', grey, 'linestyle', '--');
        plot(x1*0+1*nanmean(ILF_uf.ar.err_EC(k,:)+ ILF.decay.err_EC(k,:)), 'color', grey,'linestyle', '--');
    end
    
    ylim([-3,3]);
    
end

keyboard;
%% Generate a noise-free model

A = nanmean(p_lr.tst1(:,1)); %0.81
%B = nanmean(LR_all(:,2)); %0.38
B = nanmean(p_lr.tst1(:,2));
C = nanmean(p_lr.tst1(:,end)); %the offset

%model
xx = zeros(length(ierr_test1),1);

for qq=1:length(xx)-1
   xx(qq+1) = A*xx(qq) + B*ierr_test1(qq); %+ C; 
end

figure; hold on;
plot(xx);

plot(nanmean(res_stiff.tst1,1));

figure; plot(nanmean(norm_err.tst1,1)); %plot normalized errors

doublet_check = xx*nan;

for k=1:length(xx)-1
    doublet_check(k) = (xx(k+1) - xx(k))/norm_err.tst1(k) + (norm_err.tst1(k+1) - norm_err.tst1(k)) + (1-A)*xx(k);
    %doublet_check(k) = (xx(k+1) - xx(k))/norm_err.tst1(k) + (ierr_test1(k+1) - ierr_test1(k)) + (1-A)*xx(k);
end



%keyboard;
