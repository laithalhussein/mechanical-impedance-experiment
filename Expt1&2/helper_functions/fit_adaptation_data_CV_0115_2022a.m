function [model_data] = fit_adaptation_data_CV_0115_2022a(AR, err, mean_subtract_flag, ...
    norm_err_flag, norm_factor, gs_flag, F, nparams)

%AR -> adaptive response
%err -> errors imposed
%data are passed a matrix for each subject
%mean_subtract_flag controls whether we want to subtract the mean adaptive response
%mean_fit_flag controls whether we want to extract model parameteres based on subject-averaged data
%F is the function handle to the model we want
%nparams is the number of parameters in the model used in F (includes offset if present)
%gs_flag inidicates if we want to use grid search
%this version does leave-one-out (LOO) CV, but does so only on mean-data
%We're not so interested in determining the parameters from LOO, but rather determining the model performance

[model_data] = deal([]);
[num_trials, num_subjects] = size(AR);

%function constants
A_grid = [0.1:0.01:0.99];
B_grid = [-0.99:0.01:0];
R2_pred = NaN;
avg_pred = nan(num_trials, 1);
if nparams==2
    IG = [0.8, -0.2];
elseif nparams==4
    IG = [0.8, -0.2, 2, 0.5];
end

%demean data if desired
if mean_subtract_flag
    AR = bsxfun(@minus, AR, nanmean(AR,1));
end

AR_avg = nanmean(AR,2);
err_avg = nanmean(err,2);

AR_mat = repmat(AR_avg, [1, num_trials]);
err_mat = repmat(err_avg, [1, num_trials]);

%for LOO, I treat the sample left out as missing
AR_mat = AR_mat+diag(NaN*diag(AR_mat));
%err_mat = err_mat+diag(NaN*diag(err_mat));
%keyboard;

%before doing LOO CV, we need the test adaptation/errors that each training case will be compared to
reg_input = [err_avg, err_avg*0+1];
[btmp, ~, res] = regress(AR_avg, reg_input);
comp_test_avg = btmp(1);
AR_test_avg_cr = res/norm_factor; %impedance removed adaptation

%normalize errors if desired
test_response = AR_test_avg_cr;
if norm_err_flag
    test_input = comp_test_avg*err_avg;
else
    test_input = err_avg;
end

%perform cross-validation
for i=1:num_trials
    AR_train = AR_mat(:, i);
    err_train = err_mat(:, i);
    
    %if i==1, AR_train = AR_train(2:end); err_train = err_train(2:end); end
    
    %remove compliance for the current training set
    reg_input = [err_train, err_train*0+1];
    [btmp, ~, res] = regress(AR_train, reg_input);
    %    [~, wid] = lastwarn;
    %    warning('off', wid)
    comp_train_avg = btmp(1);
    AR_train_avg_cr = res/norm_factor;
    
    %normalize errors if desired
    tmp_response = AR_train_avg_cr;
    if norm_err_flag
        tmp_input = comp_train_avg*err_train;
    else
        tmp_input = err_train;
    end
    
    %save estimates of learning rate with a single but sensible initial guess
    p_est_avg_train_IG = nlinfit(tmp_input, tmp_response, F, IG);
    
    %optionally perform grid search
    if gs_flag
        R2gs_train_avg = nan(length(A_grid), length(B_grid));
        pgs_train_avg = cell(length(A_grid), length(B_grid));
        
        for iA = 1:length(A_grid)
            for iB = 1:length(B_grid)
                try
                    if nparams==2
                        IG_grid = [A_grid(iA), B_grid(iB)];
                    elseif nparams==4
                        IG_grid = [A_grid(iA), B_grid(iB), 2, 0.5];
                    end
                    [pgs_train_avg{iA, iB}] = nlinfit(tmp_input, tmp_response, F, IG_grid);
                    
                    %save the R2
                    model_fit_tmp = F(pgs_train_avg{iA, iB}, tmp_input);
                    R2gs_train_avg(iA,iB) = calculate_R2(tmp_response, model_fit_tmp, nparams);
                catch
                    if isempty(pgs_train_avg{iA, iB})
                        pgs_train_avg{iA, iB} = [NaN, NaN];
                        R2gs_train_avg(iA,iB) = -Inf;
                    end
                end
            end
        end
        max_R2 = max(max(R2gs_train_avg));
        [i1,i2]=find(R2gs_train_avg==max_R2);
        p_est_train_avg = pgs_train_avg{i1,i2};
        
        %impose parameter restrictions
        if nparams==2
            p_est_train_avg(1,p_est_train_avg(1)>max(A_grid)) = max(A_grid);
            p_est_train_avg(2,p_est_train_avg(2)<min(B_grid)) = min(B_grid);
        end
    end
    
    %save parameters depending on if grid search was used or not
    if gs_flag
        ptmp = p_est_train_avg;
    else
        ptmp = p_est_avg_train_IG;
    end
    
    %done estimating parameters, now we need to test them
    model_test_tmp = F(ptmp, test_input);
    
    %save the tested adapted state
    avg_pred(i) = model_test_tmp(i);    
end

% figure; hold on;
% plot(test_response, 'r');
% plot(avg_pred, 'k');

%calculate R2 and RMSE
R2_pred = calculate_R2(test_response, avg_pred, nparams);
RMSE_pred = myRMS(test_response-avg_pred);

%save the errors, the predicted adaptation, its R2, and the RMSE
model_data.pred.response = avg_pred;
model_data.pred.R2 = R2_pred;
model_data.pred.RMSE = RMSE_pred;

return