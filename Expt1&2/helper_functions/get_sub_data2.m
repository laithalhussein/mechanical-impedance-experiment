function x = get_sub_data2(dat, idx, s)
%this function returns a 3d matrix of sub x trial x sample from the dat structure
%idx can be matrix of different indices for each subject

%s is the sign of the trial

[~,num_trials] = size(idx);
[num_samples,num_subjects] = size( dat{idx(1)});

x = nan(num_subjects, num_trials,num_samples);

for i=1:num_subjects
    
    ci_sub = idx(i,:);
    
    for k=1:length(ci_sub)
         %keyboard; 
         x(i, k, :) = dat{ci_sub(k)}(:,i);
         
         if ~isempty(s)
             x(i,k,:) = x(i,k,:) * s(k);
         end

    end
end

end