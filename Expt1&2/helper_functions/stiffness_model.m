function res = stiffness_model(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + Ce(n+1)

%initialize the parameters
if isempty(params),
    A = 0.1;
    B = 0.1;
    C = 0.1;
else
    A = params(1);
    B = params(2);
    C = params(3);
end

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
     x(k+1) = A * x(k) + B * err(k) + C * err(k+1);
end

%res = x - adap;
res = x';

end