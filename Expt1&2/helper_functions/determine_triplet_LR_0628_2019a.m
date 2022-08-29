function [b,err_norm_cut,num_rem, xc] = determine_triplet_LR_0628_2019a(x, err_norm, A)

%manually check learning rate with the following:
%b = (xn - xp)/(ec / ec*) + k(en - ep) + (1-A^2)xp

%x: mean subtracted adaptation, with stiffness
%err_norm: error in units of the FF experienced

[num_subjects, nn] = size(x); %number of trials
num_rem = rem(nn,3) + 1; %since we're looking at triplets, each trial needs 2 adjacent ones

b_all = nan(num_subjects, nn - num_rem);
%b = b_all(:,2:end);

x_corrected = b_all;
x_uncorrected = b_all;

overall_LR = nan(num_subjects,1);

for i=1:num_subjects
    for j=2: nn - num_rem
        
        x_uncorrected(i,j) = (x(i,j+1) - x(i,j-1));
        stiffness_effect_tmp = -1 * (err_norm(i,j+1) - err_norm(i,j-1) );
        forgetting_effect_tmp = -(1 - (A(i)^2) )*x(i,j-1);
        
        x_corrected(i,j) = x_uncorrected(j) + stiffness_effect_tmp + forgetting_effect_tmp;
        ar_u_tmp = x_uncorrected(j)/err_norm(j);
        ar_c_tmp = x_corrected(j)/err_norm(j);
        
        b_all(i,j) = ar_c_tmp;
        
    end
    
    
    overall_LR(i) = regress(x_corrected(i,2:end)', err_norm(i,2:end-num_rem)');
    
end

b = b_all(:,2:end);
xc = x_corrected(:,2:end);

err_norm_cut = err_norm(:, [2:end-num_rem]);

% suptitle(['Trial by trial estimate of learning rate ']);


end


% nn = length(x); %number of trials
% num_rem = rem(nn,3) + 1; %since we're looking at triplets, each trial needs 2 adjacent ones
% 
% b_all = nan(1, nn - num_rem);
% b = b_all(2:end);
% 
% x_corrected = b_all;
% x_uncorrected = b_all;
% 
% figure; hold on;
% 
% for j=2: nn - num_rem
%     x_uncorrected(j) = x(j+1) - x(j-1);
%     stiffness_effect_tmp = -1 * (err_norm(j+1) - err_norm(j-1) ); %subtract stiffness effect
%     forgetting_effect_tmp = -(1 - (A^2) )*x(j-1);
%     
%     x_corrected(j) = x_uncorrected(j) + stiffness_effect_tmp + forgetting_effect_tmp;
%     ar_u_tmp = x_uncorrected(j)/err_norm(j);
%     ar_c_tmp = x_corrected(j)/err_norm(j);
%     lr_tmp = ar_c_tmp;
%     
%     b_all(j) = lr_tmp;
%     
%     %if abs(err_norm(j))==err_size, keyboard; end
%     
% end
% 
% 
% b(:) = b_all(2:end);
% 
% plot(err_norm(2:end-num_rem),b(:),'.'); %plot by normalized error
% 
% ylabel('Learning rate');
% xlabel('normalized error');
% 
% err_norm_cut = err_norm([2:end-num_rem]);
% 
% xc = x_corrected(2:end);
% xu = x_uncorrected;
% 
% suptitle(['Trial by trial estimate of learning rate ']);