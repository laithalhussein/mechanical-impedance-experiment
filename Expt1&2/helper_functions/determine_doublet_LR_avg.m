function [b,err_norm_cut,num_rem] = determine_doublet_LR_avg(x, err_norm, A)

%manually check learning rate with the following:
%b = (xn - xc)/(ec') + 1(en' - ec') + (1-A^2)xc

%this function will estimate the learning rate based on subject-averaged data

%x: mean subtracted adaptation, with stiffness (?)
%err_norm: error in units of the FF experienced

nn = length(x); %number of trials
num_rem = rem(nn,2) + 1; %since we're looking at triplets, each trial needs 2 adjacent ones

b_all = nan(1, nn - num_rem);
b = b_all(1:end);

figure; hold on;

for j=1: nn - num_rem
    
    ar_tmp = (x(j+1) - x(j)) / abs(err_norm(j));
    stiffness_effect_tmp = 1 * (err_norm(j+1) - err_norm(j) );
    forgetting_effect_tmp = (1 - A)*x(j);
    
    lr_tmp = ar_tmp+stiffness_effect_tmp+forgetting_effect_tmp;
    
    b_all(j) = lr_tmp;
    
end

b = b_all(1:end);
err_norm_cut = err_norm([1:end-num_rem]);

plot(err_norm_cut,b,'.'); %plot by normalized error

ylabel('Learning rate');
xlabel('normalized error (mm)');

suptitle(['Trial by trial estimate of learning rate ']);

end

