function [filtered_data, bad_idx] = filter_iqr_vec_1118_2021a(data, iqr_num)

[num_trials] = length(data);
bad_idx =zeros(num_trials, 1); %let 1 represent a bad trial
filtered_data = data;


mean_stat = nanmedian(data);
current_iqr = iqr(data);

if current_iqr~=0
    for nt=1:num_trials
        current_sample =  data(nt);
        
        if ~isnan(current_sample)
            if ( current_sample > mean_stat + (iqr_num * current_iqr) )  | ...
                    ( current_sample < mean_stat - (iqr_num * current_iqr) )
                
                bad_idx(nt) = 1;
                filtered_data(nt) = nan;
            end
        end
    end
end

end