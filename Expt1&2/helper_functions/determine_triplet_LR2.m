function [b,err_norm_cut,num_rem] = determine_triplet_LR2(x, err_norm, A)

%manually check learning rate with the following:
%b = (xn - xp)/(ec / ec*) + k(en - ep) + (1-A^2)xp

%x: mean subtracted adaptation, with stiffness
%err_norm: error in units of the FF experienced

nn = length(x); %number of trials
num_rem = rem(nn,3) + 1; %since we're looking at triplets, each trial needs 2 adjacent ones

b_all = nan(1, nn - num_rem);
b = b_all(2:end);

figure; hold on;

for j=2: nn - num_rem
    
    ar_tmp = (x(j+1) - x(j-1)) / abs(err_norm(j));
    stiffness_effect_tmp = -1 * (err_norm(j+1) - err_norm(j-1) );
    forgetting_effect_tmp = (1 - (A^2) )*x(j-1);
    
    lr_tmp = ar_tmp+stiffness_effect_tmp+forgetting_effect_tmp;
    
    b_all(j) = lr_tmp;
    
end

b(:) = b_all(2:end);

plot(err_norm(2:end-num_rem),b(:),'.'); %plot by normalized error

ylabel('Learning rate');
xlabel('normalized error (mm)');

err_norm_cut = err_norm([2:end-num_rem]);

suptitle(['Trial by trial estimate of learning rate ']);

end

