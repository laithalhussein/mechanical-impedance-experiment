function x = get_sub_data_0104_2022a(dat, idx)
%this function returns a 3d matrix of sub x trial x sample from the dat structure
%idx can be matrix of different indices for each subject

[num_trials, num_subjects] = size(idx);
[~,~,num_samples] = size(dat);

x = nan(num_trials, num_subjects, num_samples);
for i=1:num_subjects    
    ci_sub = idx(:,i);
    x(:,i,:) = dat(ci_sub,i,:);
end

return