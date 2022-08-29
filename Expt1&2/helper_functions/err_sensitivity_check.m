%clear all;
close all;
home;


err = sort(normrnd(0, 2.6, [1000, 1]));

%err = tand(err) * 100;

err_sens = tanh(err);


figure; hold on; 
%plot(err, err_sens, '.', 'displayname', 's = 1');

 s=[0:0.5:3];
%s = 1000;


xlim([-10,10]);


for i=1:length(s)
    
    plot(err, s(i) * tanh(err / s(i)), '-', 'displayname', ['s = ', num2str(s(i))], 'linewidth', 2);
    
end


leg1 = legend('show');
set(leg1,'Location', 'NorthWest');


xlabel('Error (deg)');
ylabel('s * tanh(e/s)');
title('Error sensitivity function');


