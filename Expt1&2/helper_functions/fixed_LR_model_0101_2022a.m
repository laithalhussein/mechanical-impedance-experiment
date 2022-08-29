function res = fixed_LR_model_0101_2022a(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + C

%x(n) is stiffness REMOVED adaptation
% C is offset

%initialize the parameters
if isempty(params),
    A = 0.8;
else
    A = params(1);
end

%restrict parameters within a specific range
A = min(max(0.1, A), 1);
%B = max(min(B, 0), -1);

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
     x(k+1) = A * x(k) - -0.1955 * err(k);
end

%res = x - adap;
res = x;

end