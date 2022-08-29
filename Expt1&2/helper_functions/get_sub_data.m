function x = get_sub_data(dat, idx)
%this function returns a 3d matrix of sub x trial x sample from the dat structure
%idx can be matrix of different indices for each subject

[~,num_trials] = size(idx);
[num_samples,num_subjects] = size( dat{idx(1)});

x = nan(num_subjects, num_trials,num_samples);

for i=1:num_subjects
    
    ci_sub = idx(i,:);
    
    for k=1:length(ci_sub)
         %keyboard; 
         x(i, k, :) = dat{ci_sub(k)}(:,i);

    end
end

end