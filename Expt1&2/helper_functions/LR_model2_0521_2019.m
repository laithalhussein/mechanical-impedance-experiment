function res = LR_model2_0521_2019(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + C

%x(n) is stiffness REMOVED adaptation
% C is offset

%initialize the parameters
if isempty(params),
    A = 0.8;
    B = -0.2;
else
    A = params(1);
    B = params(2);
end

%restrict parameters within a specific range
% A = min(max(0.1, A), 1);
% B = max(min(B, 0), -1);

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
     x(k+1) = A * x(k) - B * err(k);
end

%res = x - adap;
res = x;

return