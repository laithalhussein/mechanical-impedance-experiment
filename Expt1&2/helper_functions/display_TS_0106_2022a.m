%%% here we'll make a figure of the training schedule

ii = test2_idx(end);

tgt_tmp = tgt_all{1}([1:ii],:);

qq_ff= find(tgt_tmp(:,2)~=0);
qq_ff2 = [qq_ff,qq_ff-1,qq_ff+1];

qq_ec= find(abs(tgt_tmp(:,end))==5);
qq_ec2 = [qq_ec,qq_ec-1,qq_ec+1];

ierr_exp = tgt_all{1}([1:ii], end); %this is in degrees
ierr_exp_mm  = tand(ierr_exp) *100; %in mm

no_err_idx = find(tgt_all{1}([1:ii],7)==0);

ierr_exp_mm(no_err_idx) = NaN;

figure; plot(ierr_exp_mm, 'k.', 'markersize', 6); hold on;

plot(qq_ec2,ierr_exp_mm(qq_ec2), 'b.', 'markersize', 12);
plot(qq_ff2,ierr_exp_mm(qq_ff2), 'r.', 'markersize', 12);

%plot the EC trials as blue dots, FF trials as red dots

%mark the different phases of the experiment on the plot

yy = [-30,30];

plot(ones(1,2).*exp_epoch_all(1),yy, '--', 'color', grey);
plot(ones(1,2).*exp_epoch_all(2),yy, '--', 'color', grey);
plot(ones(1,2).*exp_epoch_all(3),yy, '--', 'color', grey);
plot(ones(1,2).*exp_epoch_all(4),yy, '--', 'color', grey);


ylabel('Imposed error (mm)');
xlabel('Trial #');

ax = gca;
ax.XTick = [30, 130, 280, 775, 1270];

xlim([0, 1290]);