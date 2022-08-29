close all;

plot(vel_FF_pa)
hold on;
%plot(vel_EC_pa)


%%
a = [0, diff(vel_FF_pa)]*200;
plot(a);

%%
qf = nanmean(squeeze(nanmedian(f_sub.FF.P,2)),1);


%% figure out coefficient

vf1 = max(vel_FF_pa) * 5;



