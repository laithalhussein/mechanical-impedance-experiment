function y = extract_adaptation_idx(x, idx)

%Here we segment the data based on the type/size of error that was experienced on the PREVIOUS trial
%If we're interested in the adaptive response on trial n based on the error of trial n-1, then we should simply find where...
%a certain error type was experienced, and then take the adaptation of the trial immediately after

[num_trials, num_subjects] = size(x);
for ksub=1:num_subjects
    
    %extract LIE data
    y.LIE.V.P3.P(:,ksub) = x(idx.LIE.V.P3.P(:,ksub)+1, ksub);
    y.LIE.V.P3.N(:,ksub) = x(idx.LIE.V.P3.N(:,ksub)+1, ksub);
    y.LIE.V.P7.P(:,ksub) = x(idx.LIE.V.P7.P(:,ksub)+1, ksub);
    y.LIE.V.P7.N(:,ksub) = x(idx.LIE.V.P7.N(:,ksub)+1, ksub);
    y.LIE.V.P0.all(:,ksub) = x(idx.LIE.V.P0.all(:,ksub)+1, ksub);
    
    %extract MIX data
    y.MIX.V.P3.P(:,ksub) = x(idx.MIX.V.P3.P(:,ksub)+1, ksub);
    y.MIX.V.P3.N(:,ksub) = x(idx.MIX.V.P3.N(:,ksub)+1, ksub);
    y.MIX.V.P7.P(:,ksub) = x(idx.MIX.V.P7.P(:,ksub)+1, ksub);
    y.MIX.V.P7.N(:,ksub) = x(idx.MIX.V.P7.N(:,ksub)+1, ksub);
    y.MIX.V.P0.all(:,ksub) = x(idx.MIX.V.P0.all(:,ksub)+1, ksub);
    
    y.MIX.C.P3.P(:,ksub) = x(idx.MIX.C.P3.P(:,ksub)+1, ksub);
    y.MIX.C.P3.N(:,ksub) = x(idx.MIX.C.P3.N(:,ksub)+1, ksub);
    y.MIX.C.P7.P(:,ksub) = x(idx.MIX.C.P7.P(:,ksub)+1, ksub);
    y.MIX.C.P7.N(:,ksub) = x(idx.MIX.C.P7.N(:,ksub)+1, ksub);
    
    %Extract HIE data
    y.HIE.C.P3.P(:,ksub) = x(idx.HIE.C.P3.P(:,ksub)+1, ksub);
    y.HIE.C.P3.N(:,ksub) = x(idx.HIE.C.P3.N(:,ksub)+1, ksub);
    y.HIE.C.P7.P(:,ksub) = x(idx.HIE.C.P7.P(:,ksub)+1, ksub);
    y.HIE.C.P7.N(:,ksub) = x(idx.HIE.C.P7.N(:,ksub)+1, ksub);
    y.HIE.V.P0.all(:,ksub) = x(idx.HIE.V.P0.all(:,ksub)+1, ksub);    
    
    %Extract all perturbation indices
    y.pert.C.all(:,ksub) = x(idx.pert.C.all(:,ksub)+1, ksub);
    y.pert.V.all(:,ksub) = x(idx.pert.V.all(:,ksub)+1, ksub);
    
end

end