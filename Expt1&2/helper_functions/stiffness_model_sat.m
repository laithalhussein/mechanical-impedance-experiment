function res = stiffness_model_sat(params, err)
%fit the training period data on the following model
%x(n+1) = Ax(n) + B * s * tanh( e(n) / s ) + C

% A is retention
% x(n) is stiffness removed adaptation
% B is the learning rate
% s is a coefficient that determines the degree of nonlinearity, the larger
% it is, the more linear the error function is, and thus more weight is given to larger errors
% e is the error
% C is the offset

%initialize the parameters
if isempty(params),
    A = 0.8;
    B = 0.07;
    s = 5; %based on simulations, this value would lead to saturation at error values close to 2.6 degrees, which is what we suspect
    C = 0.1;
else
    A = params(1);
    B = params(2);
    s = params(3);
    C = params(4);
end

N = length(err);
[x] = deal(zeros(N,1));

for k=1:N-1,
    
     x(k+1) = A * x(k) + B * s * tanh( err(k)/s ) + C;
     
     %if k==30, keyboard; end
     
end

%res = x - adap;
res = x;

end