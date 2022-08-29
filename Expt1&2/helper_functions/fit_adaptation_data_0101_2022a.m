function [p, comp, R2, model_data] = fit_adaptation_data_0101_2022a(AR, err, mean_subtract_flag, mean_fit_flag, ...
    norm_err_flag, gs_flag, F, nparams)

%AR -> adaptive response
%err -> errors imposed
%data are passed a matrix for each subject
%mean_subtract_flag controls whether we want to subtract the mean adaptive response for each subject
%mean_fit_flag controls whether we want to extract model parameteres based on subject-averaged data
%F is the function handle to the model we want
%nparams is the number of parameters in the model used in F (includes offset if present)
%gs_flag inidicates if we want to use grid search

[p, comp, R2, model_data] = deal([]);
[num_trials, num_subjects] = size(AR);

%demean data if desired
if mean_subtract_flag
    AR = bsxfun(@minus, AR, nanmean(AR,1));
end

%%%remove compliance
AR_cr = AR*NaN; %compliance removed adaptation
comp_sub = nan(num_subjects, 1); %save each subject's estimate of compliance
for k=1:num_subjects
   reg_input = [err(:,k), err(:,k)*0+1];
   [btmp, ~, res] = regress(AR(:,k), reg_input);
   [~, wid] = lastwarn;
   warning('off', wid)
   comp_sub(k) = btmp(1);
   AR_cr(:,k) = res;
end

%%%fit the model to the complaince removed adaptation using grid search
A_grid = [0.1:0.01:0.99];
B_grid = [-0.99:0.01:0];
[p_est_sub, p_est_sub_IG] = deal(nan(nparams, num_subjects));
[p_est_avg, p_est_avg_IG] = deal(NaN);
R2_sub = nan(num_subjects, 1);
R2_avg = NaN;
[sub_fit, sub_response, sub_err_norm] = deal(nan(num_trials, num_subjects));
avg_fit = nan(num_trials, 1);

for k=1:num_subjects
    
    tmp_response = AR_cr(:,k);
    if norm_err_flag
        tmp_input = comp_sub(k)*err(:,k);
    else
        tmp_input = err(:,k);
    end
    
    %save data
    sub_response(:,k) = AR_cr(:,k);
    sub_err_norm(:,k) = comp_sub(k)*err(:,k);
    
    %if the first sample of the response is NaN'd, find the first index where it isnt
    first_non_nan = 1;
    if isnan(tmp_response(1)) || isnan(tmp_input(1))
        first_non_nan = find(~isnan(tmp_response) & ~isnan(tmp_input),1,'first');
        tmp_response = tmp_response(first_non_nan:end);
        tmp_input = tmp_input(first_non_nan:end);
    end
    
    %save estimates of learning rate with a standard initial guess
    if nparams==2
        IG = [0.8, -0.2];
    elseif nparams==4
        IG = [0.8, -0.2, 2, 0.5];
    end
    p_est_sub_IG(:,k) = nlinfit(tmp_input, tmp_response, F, IG);
    
    %perform grid search
    if gs_flag
        R2gs = nan(length(A_grid), length(B_grid)); %holds R2 for each value of A and B tested
        pgs = cell(length(A_grid), length(B_grid)); %holds the corresponding parameters
        
        for iA = 1:length(A_grid)
            for iB = 1:length(B_grid)
                try
                    if nparams==2
                        IG = [A_grid(iA), B_grid(iB)];
                    elseif nparams==4
                        IG = [A_grid(iA), B_grid(iB), 2, 0.5];
                    end
                    [pgs{iA, iB}] = nlinfit(tmp_input, tmp_response, F, IG);
                    %save the R2
                    model_fit = F(pgs{iA, iB}, tmp_input);
                    R2gs(iA,iB) = calculate_R2(tmp_response, model_fit, nparams);
                catch
                    if isempty(pgs{iA, iB})
                        pgs{iA, iB} = [NaN, NaN];
                        R2gs(iA,iB) = -Inf;
                    end
                end
            end
        end
        max_R2 = max(max(R2gs));
        [i1,i2]=find(R2gs==max_R2);
        p_est_sub(:,k) = pgs{i1,i2};
        R2_sub(k) = max_R2;
        sub_fit([first_non_nan:end],k) = F(p_est_sub(:,k), tmp_input);
    end
end

%repeat if we want to estimate the model parameters based on mean data
if mean_fit_flag
    AR_avg = nanmean(AR,2);
    err_avg = nanmean(err,2);
    
%     figure;
%     plot(err_avg, AR_avg, '.');
    %keyboard;
    %remove compliance
    reg_input = [err_avg, err_avg*0+1];
    [btmp, ~, res] = regress(AR_avg, reg_input);
    %    [~, wid] = lastwarn;
    %    warning('off', wid)
    comp_avg = btmp(1);
    AR_avg_cr = res;
    
    %normalize errors if desired
    tmp_response = AR_avg_cr;
    if norm_err_flag
        tmp_input = comp_avg*err_avg;
    else
        tmp_input = err_avg;
    end
    
    %save data
    avg_response = AR_avg_cr;
    avg_err_norm = comp_avg*err_avg;
    
    %save estimates of learning rate with a standard initial guess
    if nparams==2
        IG = [0.8, -0.2];
    elseif nparams==4
        IG = [0.8, -0.2, 2, 0.5];
    end
    p_est_avg_IG = nlinfit(tmp_input, tmp_response, F, IG);
    model_fit_avg_IG = F(p_est_avg_IG, tmp_input);
    if ~gs_flag
        R2_avg = calculate_R2(tmp_response, model_fit_avg_IG, nparams);
    end
    
    %perform grid search
    if gs_flag
        R2gs_avg = nan(length(A_grid), length(B_grid));
        pgs_avg = cell(length(A_grid), length(B_grid));
        
        for iA = 1:length(A_grid)
            for iB = 1:length(B_grid)
                try
                    if nparams==2
                        IG = [A_grid(iA), B_grid(iB)];
                    elseif nparams==4
                        IG = [A_grid(iA), B_grid(iB), 2, 0.5];
                    end
                    [pgs_avg{iA, iB}] = nlinfit(tmp_input, tmp_response, F, IG);
                    
                    %save the R2
                    model_fit = F(pgs_avg{iA, iB}, tmp_input);
                    R2gs_avg(iA,iB) = calculate_R2(tmp_response, model_fit, nparams);
                catch
                    if isempty(pgs_avg{iA, iB})
                        pgs_avg{iA, iB} = [NaN, NaN];
                        R2gs_avg(iA,iB) = -Inf;
                    end
                end
            end
        end
        max_R2 = max(max(R2gs_avg));
        [i1,i2]=find(R2gs_avg==max_R2);
        p_est_avg = pgs_avg{i1,i2};
        R2_avg = max_R2;
        avg_fit = F(p_est_avg, tmp_input);        
    end
end

%impose parameter restrictions
if nparams==2
    p_est_sub(1,p_est_sub(1,:)>max(A_grid)) = max(A_grid);
    p_est_sub(2,p_est_sub(2,:)<min(B_grid)) = min(B_grid);
end

%cnstrain initial guess estimate to grid search
p_est_sub_IG(2,p_est_sub_IG(2,:)<-1) = -0.99;

p.sub.gs = p_est_sub;
p.sub.ig = p_est_sub_IG;
p.avg.gs = p_est_avg;
p.avg.ig = p_est_avg_IG;

comp.sub = comp_sub;
comp.avg = comp_avg;

R2.sub = R2_sub;
R2.avg = R2_avg;

model_data.sub.fit = sub_fit;
model_data.sub.response = sub_response;
model_data.sub.err_norm = sub_err_norm;
model_data.avg.fit = avg_fit;
model_data.avg.response = avg_response;
model_data.avg.err_norm = avg_err_norm;

return