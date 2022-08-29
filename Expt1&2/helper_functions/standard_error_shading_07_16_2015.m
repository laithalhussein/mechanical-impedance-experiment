function standard_error_shading_07_16_2015(m,s,r,N,standard_error_color)
% Plots the standard error of a time series using a shadded background
%
% m                     : array containing the mean of the data at each time point
% s                     : array contating the standard deviation of the data at each time point
% r                     : the range of the time series (ex. 1:1:100 or -50:1:50)
% N                     : Number of subjects 
% standard_error_color  : a three number array specifying the color of the shading

n=length(r);
for i=1:n-1
    X=[r(i) r(i) r(i+1) r(i+1) r(i)];
    Y=[m(i)-(s(i))/sqrt(N-1) m(i)+(s(i))/sqrt(N-1) m(i+1)+(s(i+1))/sqrt(N-1) m(i+1)-(s(i+1))/sqrt(N-1) m(i)-(s(i))/sqrt(N-1)];
    g=fill(X,Y,standard_error_color);
    h4 = findobj(g,'Type','patch');
    set(h4,'EdgeColor','none','edgealpha',0,'facecolor',standard_error_color,'facealpha',.5)
    uistack(h4, 'bottom');
    hAnnotation = get(h4,'Annotation');
    hLegendEntry = get(hAnnotation,'LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
end

return