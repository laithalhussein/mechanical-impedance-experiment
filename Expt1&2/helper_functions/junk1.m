
Nt=80;
%cfield = 'vEC';
for i=1:num_subjects, for k=1:24
        FPpre(:,i,k) = fr{idx.null.pre(i,k)}(:,i);
        FPpost(:,i,k) = fr{idx.null.post(i,k)}(:,i);
    end; end
for i=1:num_subjects, for k=1:24
        FPpre(:,i,k,2) = fr{idx.zEC.pre(i,k)}(:,i);
        FPpost(:,i,k,2) = fr{idx.zEC.post(i,k)}(:,i);
    end; end
for i=1:num_subjects, for k=1:24
        FPpre(:,i,k,3) = fr{idx.vEC.pre(i,k)}(:,i);
        FPpost(:,i,k,3) = fr{idx.vEC.post(i,k)}(:,i);
    end; end

figure; 
kk=1:1;
for k1=1:3,
for i=1:num_subjects,
    subplot(3,6,i+(k1-1)*6);hold on; 
    dFall = [FPpost(:,i,kk,k1) - FPpre(:,i,kk,k1)].*repmat(shiftdim(-FF_sign.vEC(i,kk),-1),[Nt,1,1]);
    dFmean = nanmean(dFall,3);
    dFmed = nanmedian(dFall,3);
    dFstd = nanstd(dFall,[],3);
    dFiqr = iqr(dFall,3);
    for k=kk
        
        
        
        %plot the actual pre and post forces
        
%         plot(fr{idx.null.pre(i,k)}(:,i), 'r');
%         plot(fr{idx.null.post(i,k)}(:,i), 'b');
%         plot(fr{idx.null.post(i,k)}(:,i)-fr{idx.null.pre(i,k)}(:,i), 'g','linewidth',2);
        ylim([-8 8]); 
        %plot the ideal
        plot(vel{idx.null.pre(i,k)}(:,i) * B *abs(-FF_sign.vEC(i,k)) , 'k'); hold on;
        
        %plot the AR
        dFk = dFall(:,:,k);
        plot(smooth(dFk,7), 'color', [0.5, 0, 0.5],'linewidth',0.5);
        plot(dFmean, 'color', 'g','linewidth',2);
%         plot(dFmed, 'color', 'c','linewidth',2);
%         plot(dFstd, 'color', 'r','linewidth',2);
%         plot(dFiqr/1.35, 'color', 'b','linewidth',2);
        
        grid on;
        
        %title(['AC is: ', num2str(AC.ar.null(i,k),2)] );
        
    end
    if k1==1, title(['Post-Pre: ', info.sublist{i}, ': ', info.exp_seq{i}]); end
    
end
end
