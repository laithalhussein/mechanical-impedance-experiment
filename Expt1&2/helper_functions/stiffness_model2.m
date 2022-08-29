function res = stiffness_model2(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + C

%x(n) is stiffness removed adaptation
% C is offset

%initialize the parameters
if isempty(params),
    A = 0.8;
    B = 0.07;
    C = 0.1;
else
    A = params(1);
    B = params(2);
    C = params(3);
end

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
     x(k+1) = A * x(k) + B * err(k) + C;
end

%res = x - adap;
res = x;

end