close all;
clear all;
home

% convert_imt_gz2mat2('*','C:\Users\Laith Alhussein\Desktop\Research\HSE_exp\expb\3deg_data\raw_data',...
%     'C:\Users\Laith Alhussein\Desktop\Research\HSE_exp\expb\3deg_data\mat_files',2)

%Ryan's computer
convert_imt_gz2mat2('*','C:\Users\ryanm\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\raw_data',...
    'C:\Users\ryanm\Dropbox (HNL)\Laith_files\HSE_exp\expb\3deg_data\mat_files\ALL_DATA',2)

%% crunch data
[dat,info]=run_basic_crunch_hse_shooting_1206_2021a(Cdr); %Cdr should be location of _mat folder

%save it to compile and run
%C=strcat('ALL_DATA_',info.sublist);
save(['ALL_DATA_S' num2str(size(dat.n,2)) '_' date '.mat'],'dat','info', '-v7.3');  