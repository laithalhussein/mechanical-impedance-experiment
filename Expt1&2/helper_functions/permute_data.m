function new_dat = permute_data(dat, p, n)
%p is the new order we want
%n is subject

%baseline = [1:250];
% first_condition = [251:968];
% second_condition = [969:1686];
% third_condition = [1687:2404];

exp_seq = {[251:968], [969:1686], [1687:2404]}; %note the first trial in any given block is empty
new_exp_seq = exp_seq(p);

exp_seq_tmp = exp_seq{:};
new_exp_seq_tmp = new_exp_seq{:};

dat.good(exp_seq_tmp,n) = dat.good(new_exp_seq_tmp,n);
dat.timing_data(exp_seq_tmp,n) = dat.timing_data(new_exp_seq_tmp,n);
dat.vmax(exp_seq_tmp,n) = dat.vmax(new_exp_seq_tmp,n);

dat.mtime(exp_seq_tmp,n) = dat.mtime(new_exp_seq_tmp,n);
dat.pathlen(exp_seq_tmp,n) = dat.pathlen(new_exp_seq_tmp,n);
dat.mdist(exp_seq_tmp,n) = dat.mdist(new_exp_seq_tmp,n);
dat.pathlen_rel(exp_seq_tmp,n) = dat.pathlen_rel(new_exp_seq_tmp,n);
dat.vtmax(exp_seq_tmp,n) = dat.vtmax(new_exp_seq_tmp,n);
dat.vxmax(exp_seq_tmp,n) = dat.vxmax(new_exp_seq_tmp,n);
dat.vymax(exp_seq_tmp,n) = dat.vymax(new_exp_seq_tmp,n);
dat.ra(exp_seq_tmp,n) = dat.ra(new_exp_seq_tmp,n);
dat.tvtmax(exp_seq_tmp,n) = dat.tvtmax(new_exp_seq_tmp,n);
dat.n(exp_seq_tmp,n) = dat.n(new_exp_seq_tmp,n);

num_samples = size(dat.vxr{2}(:,1),1);
num_trials = size(dat.vxr,2)

%pad all empty trials with nans

for k=1:num_trials
    
   if isempty(dat.vxr)
       
   end
    
    
    
    
    
end


for q=1:length(exp_seq_tmp)
    
    if ~isempty(dat.vxr{exp_seq_tmp(q)})
        dat.vxr{exp_seq_tmp(q)}(:,n) = dat.vxr{new_exp_seq_tmp(q)}(:,n);
        dat.vyr{exp_seq_tmp(q)}(:,n) = dat.vyr{new_exp_seq_tmp(q)}(:,n);
        dat.pxr{exp_seq_tmp(q)}(:,n) = dat.pxr{new_exp_seq_tmp(q)}(:,n);
        dat.pyr{exp_seq_tmp(q)}(:,n) = dat.pyr{new_exp_seq_tmp(q)}(:,n);
    end
    
    if isempty(dat.fyr{exp_seq_tmp(q)})
        
        
        
        dat.fyr{exp_seq_tmp(q)}(:,n) = dat.fyr{new_exp_seq_tmp(q)}(:,n);
    end

end

new_dat = dat;


end