function [AC, rsq] = calculate_AC(F_data, B, V_data)
%data comes in as a 3D matrix: (subject, trial, sample)

[nsubject,ntrial, nsample] = size(F_data);
AC = nan(nsubject,ntrial);
rsq = nan(nsubject,ntrial);


% if ntrial>ns
%     F_data = F_data';
%     V_data = V_data';
%     [ns,ntrial] = size(F_data);
% end


for k=1:nsubject
    for i=1:ntrial
        
        F = squeeze(F_data(k,i,:)) * sign(B(i));
        V = squeeze(V_data(k,i,:)) * sign(B(i));
        %     ideal_F = V * abs(B); %use if we want all the same signs for learning curves
        ideal_F = V * (abs(B(i)));
        
        %plot(ideal_F); hold on; plot(F);
        %keyboard;
        if ~(sum(isnan(F)) == nsample)
            reg_input = [ideal_F*0+1, ideal_F];
            
            [b,~,~,~,s] = regress(F, reg_input);
            
            AC(k,i) = b(2);
            rsq(k,i) = s(1);
            
        end
        
    end
end

 %keyboard;

end