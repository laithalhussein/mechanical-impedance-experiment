function [ h, orient_deg ] = plot_error_ellipse(x,y,col,LW, show_CI)
%need to plot an ellipse where the axes are defined by the noise we have in the data

%how many errors do we want each axis of the ellipse to signify
num_errors = 1;
if strcmp(show_CI, 'show_CI')
    num_errors = 1/sqrt(length(x))*1.96;
end

c = nancov([x(:),y(:)]);

%get the eigenvalues/vectors
[evec,eval] = eig(c);

%start to define the ellipse points
a = [0:360]' * pi/180;
e_points = [cos(a),sin(a)]*sqrt(eval)*evec'*num_errors;

h = plot(e_points(:,1)+nanmean(x), e_points(:,2)+nanmean(y),'color',col,'linewidth',LW); %

%h = plot(exy(:,1)+x_center, exy(:,2)+y_center);

uistack(h, 'bottom');
hg=get(h,'Annotation');
hLegendEntry = get(hg,'LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off');

% find the orientation
[max_eval,NDX]=max(max(abs(eval)));

orient_deg = 180/pi*atan2(evec(2,NDX),evec(1,NDX));

% get the area???
H=sqrt(eval)*evec'*num_errors;
r_1=sqrt(H(1,1)^2+H(1,2)^2);
r_2=sqrt(H(2,1)^2+H(2,2)^2);
R=[r_1,r_2];

% ellipse_area = r_1*r_2*pi;

hold on

%we could show the eigenvector itself...
%h1=plot([mean(x),mean(x)+R(NDX)*cosd(orient_deg)],[mean(y),mean(y)+R(NDX)*sind(orient_deg)],'LineWidth',3,'color','g');

shg
return

