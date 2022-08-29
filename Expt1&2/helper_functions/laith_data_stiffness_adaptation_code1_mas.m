function laith_data_stiffness_adaptation_code1_mas(v2,p2)

 %make_plot3(v2.f_sub,'E3 Fx vFF:  ')

%keyboard

make_plot1(v2.f_sub,'E3 Fx vFF:  ', 'X-Force (N)')
%make_plot1(p2.f_sub,'E4 Fx pFF:  ', 'X-Force (N)')

make_plot2(v2.ypos_sub, v2.xpos_sub,'E3 x-y vFF:  ', 'Y-Position (mm)','X-Position (mm)')
%make_plot2(p2.ypos_sub, p2.xpos_sub,'E4 x-y pFF:  ', 'Y-Position (mm)','X-Position (mm)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make_plot2(v2.f_sub, v2.ypos_sub, 'E3 Fx vFF:  ', 'X-Force (N)','Y-Position (mm)')
% make_plot2(p2.f_sub, p2.ypos_sub, 'E4 Fx pFF:  ', 'X-Force (N)','Y-Position (mm)')

% make_plot2(v2.xpos_sub, v2.ypos_sub,'E3 y-x vFF:  ', 'X-Position (mm)','Y-Position (mm)')
% make_plot2(p2.xpos_sub, p2.ypos_sub,'E4 y-x pFF:  ', 'X-Position (mm)','Y-Position (mm)')

% make_plot1(v2.xpos_sub,'E3 x-pos vFF:  ', 'X-Position (mm)')
% make_plot1(p2.xpos_sub,'E4 x-pos pFF:  ', 'X-Position (mm)')
% 
% make_plot1(v2.ypos_sub,'E3 y-pos vFF:  ', 'Y-Position (mm)')
% make_plot1(p2.ypos_sub,'E4 y-pos pFF:  ', 'Y-Position (mm)')
end

function make_plot1(q,label1,label2)
purple = [147,112,219]/255;
ep = squeeze(nanmedian(q.EC.P,2))';
en = squeeze(nanmedian(q.EC.N,2))';
fp = squeeze(nanmedian(q.FF.N,2))';
fn = squeeze(nanmedian(q.FF.P,2))';
mep = mean(ep,2); %subject avg
men = mean(en,2);
mfp = mean(fp,2); 
mfn = mean(fn,2);
%figure;

%plot in terms of time instead of sample
N = size(ep,2);
tt = [0:length(mfn)-1]*0.005;

figure; hold on; grid; title([label1,'+/- perturbations (thin/thick = EC/FF)']);  
plot(tt,ep,'b-'); plot(tt,en,'r--'); plot(tt,fp,'b-','linewidth',2); plot(tt,fn,'r--','linewidth',2);
ylim([-20,20]);

figure; hold on; grid; title('+/- perturbations');
plot(tt,mep,'b-'); plot(tt,men,'r--'); plot(tt,mfp,'linewidth',2); plot(tt,mfn,'r--','linewidth',2);

%add the standard error
standard_error_shading_07_16_2015(mep,std(ep,[],2),tt,N,'b');
standard_error_shading_07_16_2015(men,std(en,[],2),tt,N,'r');
standard_error_shading_07_16_2015(mfp,std(fp,[],2),tt,N,'b');
standard_error_shading_07_16_2015(mfn,std(fn,[],2),tt,N,'r');

ylim([-10,10]);

figure; hold on; grid; title('+/- perturbation diff'); 
plot(tt,ep-en,'-', 'color', purple); plot(tt,fp-fn,'-', 'color', purple,'linewidth',2);
ylim([-10,10]);

figure; hold on; grid; title('+/- perturbation diff'); 
plot(tt,mep-men,'-', 'color', purple); plot(tt,mfp-mfn,'-', 'color', purple,'linewidth',2);
ylim([-2,10]);

standard_error_shading_07_16_2015(mep-men,std(ep-en,[],2),tt,N,purple);
standard_error_shading_07_16_2015(mfp-mfn,std(fp-fn,[],2),tt,N,purple);

%subplot(221); ylabel(label2); xlabel('Time (ms)');
end

function make_plot2(q,q2,label1,label2,label3)
purple = [147,112,219]/255;
ep = squeeze(nanmedian(q.EC.P,2))';
en = squeeze(nanmedian(q.EC.N,2))';
fp = squeeze(nanmedian(q.FF.N,2))';
fn = squeeze(nanmedian(q.FF.P,2))';
mep = mean(ep,2); %subject avg
men = mean(en,2);
mfp = mean(fp,2); 
mfn = mean(fn,2);

ep2 = squeeze(nanmedian(q2.EC.P,2))';
en2 = squeeze(nanmedian(q2.EC.N,2))';
fp2 = squeeze(nanmedian(q2.FF.N,2))';
fn2 = squeeze(nanmedian(q2.FF.P,2))';
mep2 = mean(ep2,2); %subject avg
men2 = mean(en2,2);
mfp2 = mean(fp2,2); 
mfn2 = mean(fn2,2);

N = size(ep,2);

figure; 
hold on; grid; title([label1,'+/- perturbations (thin/thick = EC/FF)']);  
plot(ep2,ep,'b-'); plot(en2,en,'r--'); plot(fp2,fp,'b-','linewidth',2); plot(fn2,fn,'r--','linewidth',2); 
xlim([-45,45]); ylim([-25,150]);


figure; hold on; grid; title('+/- perturbations');     
plot(mep2,mep,'b-'); plot(men2,men,'r--'); plot(mfp2,mfp,'linewidth',2); plot(mfn2,mfn,'r--','linewidth',2);
xlim([-30,30]); ylim([-25,150]);

standard_error_shading_07_16_2015(mep,std(ep,[],2),mep2,N,'b');
standard_error_shading_07_16_2015(men,std(en,[],2),men2,N,'r');
standard_error_shading_07_16_2015(mfp,std(fp2,[],2),mfp2,N,'b');
standard_error_shading_07_16_2015(mfn,std(fn2,[],2),mfn2,N,'r');

%if x is on the y-axis, use the following, otherwise see below
% subplot(222); hold on; grid; title('+/- perturbation diff'); plot(ep2,ep-en,'-', 'color', purple); plot(fp2,fp-fn,'-', 'color', purple,'linewidth',2);
% subplot(224); hold on; grid; title('+/- perturbation diff'); plot(mep2,mep-men,'-', 'color', purple); plot(mfp2,mfp-mfn,'-', 'color', purple,'linewidth',2);

figure; hold on; grid; title('+/- perturbation diff'); 
plot(ep2-en2,ep,'-', 'color', purple); plot(fp2-fn2,fp,'-', 'color', purple,'linewidth',2);
xlim([-45,45]); ylim([-25,150]);


figure; hold on; grid; title('+/- perturbation diff'); 
plot(mep2-men2,mep,'-', 'color', purple); plot(mfp2-mfn2,mfp,'-', 'color', purple,'linewidth',2);
xlim([-30,30]); ylim([-25,150]);

standard_error_shading_07_16_2015(mep,std(ep2-en2,[],2),mep2-men2,N,purple);
standard_error_shading_07_16_2015(mfp,std(fp2-fn2,[],2),mfp2-mfn2,N,purple);

%subplot(221); ylabel(label2); xlabel(label3);
end



function make_plot3(q,label1)
purple = [147,112,219]/255;
%sN = 10;
ep = squeeze(nanmean(q.EC.P,2))';
en = squeeze(nanmean(q.EC.N,2))';
fp = squeeze(nanmean(q.FF.N,2))';
fn = squeeze(nanmean(q.FF.P,2))';
epa = squeeze(nanmean(q.ECpost.P-q.ECpre.P,2))';
ena = squeeze(nanmean(q.ECpost.N-q.ECpre.N,2))';
fpa = squeeze(nanmean(q.FFpost.N-q.FFpre.N,2))';
fna = squeeze(nanmean(q.FFpost.P-q.FFpre.P,2))';

sN = sqrt(size(ep,2));

tt = [0:size(ep,1)-1]*0.005;


%figure;
figure; hold on; title([label1, 'EC perturbation and adaptive response']); plot(tt,ep-en,'b'); plot(tt,-(epa-ena),'r');
ylim([-10,10]);

figure; hold on; title('FF perturbation and adaptive response'); plot(tt,fp-fn,'b'); plot(tt,-(fpa-fna),'r');
ylim([-10,10]);

figure;hold on; plot(tt,mean(ep-en,2),'b'); plot(tt,mean((epa-ena),2),'r');

standard_error_shading_07_16_2015(mean(ep-en,2),std(ep-en,[],2),tt,sN,'b');
standard_error_shading_07_16_2015(mean(epa-ena,2),std(epa-ena,[],2),tt,sN,'r');
ylim([-2,10]);

figure; hold on; plot(tt,mean(fp-fn,2),'b'); plot(tt,mean((fpa-fna),2),'r');

standard_error_shading_07_16_2015(mean(fp-fn,2),std(fp-fn,[],2),tt,sN,'b');
standard_error_shading_07_16_2015(mean(fpa-fna,2),std(fpa-fna,[],2),tt,sN,'r');
ylim([-2,10]);

% subplot(223); hold on; plot(median(ep-en,2),':b'); plot(median((epa-ena),2),':r');
% subplot(224); hold on; plot(median(fp-fn,2),':b'); plot(median((fpa-fna),2),':r');

% subplot(223); grid;hold on; errorbar(mean(ep-en,2),std(ep-en,[],2)/sN, 'b'); errorbar(mean((epa-ena),2), std((epa-ena),[],2)/sN,'r');
% subplot(224); grid;hold on; errorbar(mean(fp-fn,2),std(fp-fn,[],2)/sN, 'b'); errorbar(mean((fpa-fna),2), std((fpa-fna),[],2)/sN, 'r');
% subplot(223); hold on; plot(median(ep-en,2),':b'); plot(median((epa-ena),2),':r');
% subplot(224); hold on; plot(median(fp-fn,2),':b'); plot(median((fpa-fna),2),':r');
%keyboard;
%subplot(221);  xlabel('Time step (dt=5ms)');
end



% figure; plot(squeeze(mean(v2.xpos_sub.EC.P,2))'); hold on; plot(squeeze(mean(v2.xpos_sub.EC.N,2))','-'); shg
% figure; plot(squeeze(mean(v2.xpos_sub.FF.P,2))'); hold on; plot(squeeze(mean(v2.xpos_sub.FF.N,2))','--'); shg
% 
% purple = [147,112,219]/255;
% ep = squeeze(nanmedian(v2.xpos_sub.EC.P,2))';
% vxen = squeeze(nanmedian(v2.xpos_sub.EC.N,2))';
% vxfp = squeeze(nanmedian(v2.xpos_sub.FF.P,2))';
% vxfn = squeeze(nanmedian(v2.xpos_sub.FF.N,2))';
% mvxep = mean(vxep,2); %subject avg
% mvxen = mean(vxen,2);
% mvxfp = mean(vxfp,2); 
% mvxfn = mean(vxfn,2);
% figure; 
% subplot(221); hold on; grid; title('E3: +/- perturbations (thin/thick = EC/FF)');     plot(vxep,'b-'); plot(vxen,'r--'); plot(vxfp,'b-','linewidth',2); plot(vxfn,'r--','linewidth',2); 
% subplot(223); hold on; grid; title('+/- perturbations');     plot(mvxep,'b-'); plot(mvxen,'r--'); plot(mvxfp,'linewidth',2); plot(mvxfn,'r--','linewidth',2);
% subplot(222); hold on; grid; title('+/- perturbation diff'); plot(vxep-vxen,'-', 'color', purple); plot(vxfp-vxfn,'-', 'color', purple,'linewidth',2);
% subplot(224); hold on; grid; title('+/- perturbation diff'); plot(mvxep-mvxen,'-', 'color', purple); plot(mvxfp-mvxfn,'-', 'color', purple,'linewidth',2);
% subplot(221); ylabel('x-displacement (mm)'); xlabel('Time step (dt=5ms)');
% 
% vxep = squeeze(nanmedian(p2.xpos_sub.EC.P,2))';
% vxen = squeeze(nanmedian(p2.xpos_sub.EC.N,2))';
% vxfp = squeeze(nanmedian(p2.xpos_sub.FF.P,2))';
% vxfn = squeeze(nanmedian(p2.xpos_sub.FF.N,2))';
% mvxep = mean(vxep,2); %subject avg
% mvxen = mean(vxen,2);
% mvxfp = mean(vxfp,2); 
% mvxfn = mean(vxfn,2);
% figure; 
% subplot(221); hold on; grid; title('E4: +/- perturbations (thin/thick = EC/FF)');     plot(vxep,'b-'); plot(vxen,'r--'); plot(vxfp,'b-','linewidth',2); plot(vxfn,'r--','linewidth',2); 
% subplot(223); hold on; grid; title('+/- perturbations');     plot(mvxep,'b-'); plot(mvxen,'r--'); plot(mvxfp,'linewidth',2); plot(mvxfn,'r--','linewidth',2);
% subplot(222); hold on; grid; title('+/- perturbation diff'); plot(vxep-vxen,'-', 'color', purple); plot(vxfp-vxfn,'-', 'color', purple,'linewidth',2);
% subplot(224); hold on; grid; title('+/- perturbation diff'); plot(mvxep-mvxen,'-', 'color', purple); plot(mvxfp-mvxfn,'-', 'color', purple,'linewidth',2);
% subplot(221); ylabel('x-displacement (mm)'); xlabel('Time step (dt=5ms)');
% 
% 
% figure; subplot(111); hold on; 
% plot(vxep); plot(vxen,'--'); plot(vxfp,'linewidth',2); plot(vxfn,'--','linewidth',2);
% 
% figure; subplot(111); hold on; 
% plot(vxep-vxen); plot(vxfp-vxfn,'linewidth',2);
% 
% 
% figure; plot(mean(squeeze(nanmedian(v2.xpos_sub.EC.P,2))',2)); hold on; plot(mean(squeeze(nanmedian(v2.xpos_sub.EC.N,2))',2),'--'); shg
% hold on; plot(mean(squeeze(nanmedian(v2.xpos_sub.FF.P,2))',2), 'linewidth',2); hold on; plot(mean(squeeze(nanmedian(v2.xpos_sub.FF.N,2))',2),'--', 'linewidth',2); shg
% 
% figure; 
% subplot(221); hold on; grid; title('+/- perturbations');     plot(vxep); plot(vxen,'--'); plot(vxfp,'linewidth',2); plot(vxfn,'--','linewidth',2); 
% subplot(222); hold on; grid; title('+/- perturbation diff'); plot(vxep-vxen); plot(vxfp-vxfn,'linewidth',2);
% subplot(223); hold on; grid; title('+/- perturbations');     plot(mvxep); plot(mvxen,'--'); plot(mvxfp,'linewidth',2); plot(mvxfn,'--','linewidth',2);
% subplot(224); hold on; grid; title('+/- perturbation diff'); plot(mvxep-mvxen); plot(mvxfp-mvxfn,'linewidth',2);
