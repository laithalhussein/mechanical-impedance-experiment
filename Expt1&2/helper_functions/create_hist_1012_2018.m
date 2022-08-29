function create_hist_1012_2018(data,color,b)

if isempty(b)
    bins = 18;
else 
    bins = b;
end

alpha_val = 0.6;

if isempty(data), return, end

%hh=histogram(data,bins,'facecolor',color,'facealpha',alpha_val, 'edgecolor','none');
hh=histogram(data,bins, 'DisplayStyle', 'stairs', 'linewidth', 1);

hh.EdgeColor = color;

%add vertical line @ mean
mu = nanmean(data);
%mu = nanmedian(data);
YL = ylim;

h1=plot([mu, mu], [YL(1), YL(2)], 'color',color, 'linestyle','--','linewidth',3,'Displayname',num2str(mu));
% leg1=legend(h1);
% set(leg1,'Location','Best');

% plot([T, T], [YL(1), YL(2)], 'color',ones(1,3) * 0, 'linestyle','-','linewidth',3); %meant for median
% leg1 = legend([num2str(mu)]);
% set(leg1,'Location','Best');

%axis tight;
end