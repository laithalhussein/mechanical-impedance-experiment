function res = stiffness_model_refresh(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + D(n)

%x(n) is stiffness removed adaptation
%D(n) is a vector that should include offsets for when there are discontinuities in the data

%initialize the parameters

A = params(1);
B = params(2);
%C = params(3);
D = params(3:end);

% err = xi{1};
% idx_disc = xi{2};

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
     x(k+1) = A * x(k) + B * err(k) + D(k);
end

%res = x - adap;
res = x;

end