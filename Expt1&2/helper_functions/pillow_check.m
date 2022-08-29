clear all;
close all;
home;

cd('C:\Users\Laith Alhussein\Desktop\Research\HSE_exp\data\no_washout\shooting_data\mat_files\ALL_data');

Cdr=cd; %Make sure we're in the folder containing the grouped data

hse_filename=ls('ALL_DATA_*'); %SM=shooting movements...
hse_data=load(strcat(Cdr,'\',hse_filename));
dat=hse_data.dat;
info=hse_data.info;


































