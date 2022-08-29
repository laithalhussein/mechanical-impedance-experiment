
close all;
home;

win = [-79:0];

%check subject 1, (ZLV), 0.7796 null condition
%is the 23 pre trial = 831 for this participant? -> yes...
%this corresponds to 1662 in original tgt file

sub_id = 4;
trial_id = 832; %

tgt_all_tmp = info.tgt_dat(sub_id,:);
tgt_all = vertcat(tgt_all_tmp{:});

%FF_sign = tgt_all(trial_id-2,2);

f_tmp1 = dat.fyr{trial_id}(:,sub_id);
qv = dat.vxr{trial_id}(:,sub_id);
qp = dat.pxr{trial_id}(:,sub_id);

pillow_idx = find( qp > 0.1, 1, 'first') - 1;


p = qp(pillow_idx + win);
v = qv(pillow_idx + win);
F1 = f_tmp1(pillow_idx + win);



figure; hold on; plot(v); plot(p);

idealF = v * 7.5; %* FF_sign;
figure; hold on; plot(idealF,'k');

plot(F1,'b'); 

b = regress(F1, [idealF, idealF*0+1])




% f_tmp2 = dat.fyr{trial_id-4}(:,sub_id);
% F2 = f_tmp2(pillow_idx + win);
% plot(F2,'r');