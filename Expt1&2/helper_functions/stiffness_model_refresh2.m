function res = stiffness_model_refresh2(params, xi)
%fit the training period data on the following model
%x(n+1) = Ax(n) + Be(n) + D(n)

%x(n) is stiffness removed adaptation
%D(n) is a vector that should include offsets for when there are discontinuities in the data

%initialize the parameters

A = params(1);
B = params(2);
%C = params(3);
D = params(3:end);

err = xi{1};
idx_disc = xi{2};

N = length(err);
[x] = deal(zeros(N,1));

c = 1;
for k=1:N-1,
    
    try

        if ismember(k,idx_disc)
            x(k+1) = A * x(k) + B * err(k) + D(c);
            c = c+1;
        else
            x(k+1) = A * x(k) + B * err(k);
        end
        
    catch
        keyboard;
    end
    
end

%res = x - adap;
res = x;

end