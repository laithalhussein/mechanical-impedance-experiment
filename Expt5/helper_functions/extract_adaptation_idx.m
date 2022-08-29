function y = extract_adaptation_idx(x, idx)

%Here we segment the data based on the type/size of error that was experienced on the PREVIOUS trial
%If we're interested in the adaptive response on trial n based on the error of trial n-1, then we should simply find where...
%a certain error type was experienced, and then take the adaptation of the trial immediately after

[num_trials, num_subjects] = size(x);
for ksub=1:num_subjects
    
    %extract NH data
    y.NH.V.P0.all(:,ksub) = x(idx.NH.V.P0.all(:,ksub)+1, ksub);
    y.NH.V.P3.P(:,ksub) = x(idx.NH.V.P3.P(:,ksub)+1, ksub);
    y.NH.V.P3.N(:,ksub) = x(idx.NH.V.P3.N(:,ksub)+1, ksub);
    y.NH.V.P7.P(:,ksub) = x(idx.NH.V.P7.P(:,ksub)+1, ksub);
    y.NH.V.P7.N(:,ksub) = x(idx.NH.V.P7.N(:,ksub)+1, ksub);
    
    y.NH.C.P3.P(:,ksub) = x(idx.NH.C.P3.P(:,ksub)+1, ksub);
    y.NH.C.P3.N(:,ksub) = x(idx.NH.C.P3.N(:,ksub)+1, ksub);
    y.NH.C.P7.P(:,ksub) = x(idx.NH.C.P7.P(:,ksub)+1, ksub);
    y.NH.C.P7.N(:,ksub) = x(idx.NH.C.P7.N(:,ksub)+1, ksub);
    
    %extract HI data
    y.HI.P3.P(:,ksub) = x(idx.HI.P3.P(:,ksub)+1, ksub);
    y.HI.P3.N(:,ksub) = x(idx.HI.P3.N(:,ksub)+1, ksub);
    y.HI.P7.P(:,ksub) = x(idx.HI.P7.P(:,ksub)+1, ksub);
    y.HI.P7.N(:,ksub) = x(idx.HI.P7.N(:,ksub)+1, ksub);
    
    %Extract HIE data
    y.HD.P3.P(:,ksub) = x(idx.HD.P3.P(:,ksub)+1, ksub);
    y.HD.P3.N(:,ksub) = x(idx.HD.P3.N(:,ksub)+1, ksub);
    y.HD.P7.P(:,ksub) = x(idx.HD.P7.P(:,ksub)+1, ksub);
    y.HD.P7.N(:,ksub) = x(idx.HD.P7.N(:,ksub)+1, ksub);
    
    %Extract all perturbation indices
    y.pert.C.all(:,ksub) = x(idx.pert.C.all(:,ksub)+1, ksub);
    y.pert.V.all(:,ksub) = x(idx.pert.V.all(:,ksub)+1, ksub);
    
end

end