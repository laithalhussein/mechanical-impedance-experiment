function [b,err_norm_cut,num_rem,xc,xu] = determine_triplet_LR_avg_ns(x, err_norm, A, esize, triplet_check, triplet_check2)

%manually check learning rate with the following:
%b = (xn - xp)/(ec / ec*) + k(en - ep) + (1-A^2)xp

%this function will estimate the learning rate based on subject-averaged data

%x: mean subtracted adaptation, with stiffness
%err_norm: error in units of the FF experienced

nn = length(x); %number of trials
num_rem = rem(nn,3) + 1; %since we're looking at triplets, each trial needs 2 adjacent ones

b_all = nan(1, nn - num_rem);
b = b_all(2:end);

x_corrected = b_all;
x_uncorrected = b_all;

figure; hold on;

for j=2: nn - num_rem
    x_uncorrected(j) = x(j+1) - x(j-1);
    stiffness_effect_tmp = 0; %subtract stiffness effect
    forgetting_effect_tmp = -(1 - (A^2) )*x(j-1);
    
    x_corrected(j) = x_uncorrected(j) + stiffness_effect_tmp + forgetting_effect_tmp;
    ar_u_tmp = x_uncorrected(j)/err_norm(j); %get rid of absolute value
    ar_c_tmp = x_corrected(j)/err_norm(j); %get rid of absolute value
    lr_tmp = ar_c_tmp; 
    
    b_all(j) = lr_tmp;
    
    
    if j==4 %err_norm(j) == esize
        
       % keyboard;
       
    end
    
end

b(:) = b_all(2:end);

plot(err_norm(2:end-num_rem),b(:),'.'); %plot by normalized error

ylabel('Learning rate');
xlabel('normalized error (mm)');

err_norm_cut = err_norm([2:end-num_rem]);

xc = x_corrected(2:end);
xu = x_uncorrected;

suptitle(['Trial by trial estimate of learning rate ']);

end

