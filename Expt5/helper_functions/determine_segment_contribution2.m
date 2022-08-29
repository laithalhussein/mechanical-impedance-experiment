function bsub = determine_segment_contribution2(v, c, hd, hi)
%we're expecting arounf 17% vs 83% for 1st vs second half effect
s = 1;
nsub = length(v);

%containate the data for design and response matrix
y_sub = [v; c; hd; hi];
%encode the presence of the VC in first half vs second half (first columnis offset
X_tmp = [1 0 0;...
    1, 1, 1;...
    1, 1, 0;...
    1, 0, 1];

%regress based on the mean first, to make sure we have this making sense
b_avg = regress(mean(y_sub,2), X_tmp);

%try for every subject
bsub = nan(nsub,3);
for k=1:nsub
   bsub(k,:) = regress(y_sub(:,k), X_tmp);
   %bsub(k) = min(max(bsub(k),0),1);
end
%mean(bsub,1)
return