%%% here we'll make a figure of the training schedule

tgt_tmp = tgt_all{1};
ierr_exp = tgt_tmp(:, end); %this is in degrees
no_err_idx = find(tgt_tmp(:,7)==0);
ierr_exp(no_err_idx) = NaN;

ierr_test1 = ierr_exp(test1_idx);
ierr_test2 = ierr_exp(test2_idx);

ff_pert= find(tgt_tmp(test2_idx,2)~=0);
ff_trip = [ff_pert-1;ff_pert+1];
ff_size = tgt_tmp(test2_idx(ff_pert),2)*B;

fc_pert= find(abs(tgt_tmp(test1_idx,end))==5);
fc_trip = [fc_pert;fc_pert-1;fc_pert+1];

trials_test = [1:length(ierr_test1)];

%make 2 separate plots: (1) FC test period, and (2) FF test period

%% FC test
figure; hold on;
plot(trials_test, ierr_test1, 'color', 'k', 'markersize', 6);
plot(fc_trip,ierr_test1(fc_trip), 'color', 'b', 'linestyle', 'none', 'marker', '.', 'markersize', 12);

ylabel('Imposed FC error (deg)');
xlabel('Trial number');

ax = gca;
ax.XTick = [0, 500];
ax.YTick = [-8, 0, 8];

%% FF test
figure; hold on;
plot(trials_test, ierr_test2, 'color', 'k', 'markersize', 6);
plot(ff_trip,ierr_test2(ff_trip), 'color', 'r', 'linestyle', '-', 'marker', '.', 'markersize', 12);
%plot(ff_pert, ff_size*0, 'color', 'r', 'linestyle', 'none', 'marker', 'p', 'markersize', 12);
ylabel('Imposed FC error (deg)');
ax = gca;
ax.YTick = [-8, 0, 8];
ax.XTick = [0, 500];
xlabel('Trial number');



% yyaxis left
% plot(trials_test, ierr_test2, 'color', 'k', 'markersize', 6);
% plot(ff_trip,ierr_test2(ff_trip), 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 12);
% ylabel('Imposed error (deg)');
% ax = gca;
% ax.YTick = [-8, 0, 8];
% 
% yyaxis right
% ylabel('FF strength (Ns/m)');
% plot(ff_pert,ff_size, 'color', 'r', 'linestyle', 'none', 'marker', '.', 'markersize', 12);
% ylim([-B-1, B+1])
% ax = gca;
% ax.YTick = [-B, 0, B];
% 
% xlabel('Trial number');
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
% ax.XTick = [0, 500];
