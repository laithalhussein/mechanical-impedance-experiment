function [ R2,p ] = calculate_R2(original_data, model_data,dof )
%calculates the R^2 and optionally returns the significance
% dof (degrees of freedom) INCLUDES a constant offset, i.e. y=mx+b has 2 
 
if nargin<3
    sig_flag = 0;
else
    sig_flag = 1;
end
 
good_idx = ~isnan(original_data) & ~isnan(model_data);
original_data = original_data(good_idx);
model_data = model_data(good_idx);
 
res = original_data - model_data;
SS_r = sum(res.^2);
SS_total = var(original_data)*(length(original_data)-1);
 
R2 = 1-SS_r/SS_total;
 
if(sig_flag)
    p_tmp = 1;
    n = length(original_data);
    F = (n-dof)/(dof-p_tmp)*(1/(1-R2)-1);
    p = 1-fcdf(F,dof-p_tmp,n-dof);
else
    p = nan;
end
 
return