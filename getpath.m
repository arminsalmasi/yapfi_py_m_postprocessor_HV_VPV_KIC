%% clear work environment 
clear variables;
close all;
clc;
current_path= pwd ;  work_path = uigetdir(current_path);
clear current_path

%% get the wok_path
%work_path = uigetdir(pwd);
%work_path = 'C:\Users\armin\Documents\Mycodes\YAPFI-MATLAB-HardNess-2018\test2D';

%% next script

disp('!!!!!!!!!!!!<<get path end>>!!!!!!!!!!!!')
disp(work_path)
readprofiles

