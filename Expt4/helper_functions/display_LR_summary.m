function h = display_LR_summary(y, err, color_order, x_labels, y_label, fig_title, ymm)
%summarizes the results with bar plots

h = figure; hold on;
xb_spacing = 0.55;

num_conditions = length(y);

xb1 = linspace(1,xb_spacing*num_conditions+xb_spacing,num_conditions);

for k=1:num_conditions
%     barwitherr(err(k), xb1(k), -1*imp_LR_avg(k), xb_spacing-0.05, 'Edgecolor', 'k',...
%         'FaceColor', color_order(k,:),'Linewidth', 2); hold on;
    barwitherr(NaN, xb1(k), y(k), xb_spacing-0.05, 'Edgecolor', 'k',...
        'FaceColor', color_order(k,:),'Linewidth', 2); hold on;
    display_xy_error(xb1(k), y(k), [], err(k), 'color', 'k');
end

%Optionally plot the data on top of the bar graphs as well
% for ks=1:length(session_flds)
%     ef = env_flds_{ks};
%     plot_data_on_bar(xb1(ks), LR_reg_sub_all.(ef)(:,1), 'Marker', 'o', 'Color', 'k', 'Linewidth',1.5); hold on;
% end

xticks([xb1(:)]);
xticklabels(x_labels);
xtickangle(45)

set(gca,'layer','top','box','off');
ylabel(y_label);
title(fig_title);

if isempty(ymm)
    ylim([0,1]);    
else
    ylim(ymm);
end
%grid on;
ax = gca;

ylim([-0.05,0.6]);
%grid on;
ax.YTick = [0, 0.2, 0.4, 0.6];

return