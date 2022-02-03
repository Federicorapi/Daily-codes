clear all
clc
close all
%% 
[File, Path] = uigetfile('*.h5');   %% Insert file name
%I_Ml_on_SiO2_0.3xname = 'PL_PAS_FIN0.7V_resolution_0.005V_0.4s_exposure_spectro750_excitation632_20µW_Olympus_002.h5_001.h5';      %% Insert file name
name = strcat(Path,File);
wvl = h5read(name,'/WL');
%ev=1240/wvl;
data = h5read(name,'/Data');
spectra= data;


figure()
plot(wvl,spectra)
xlabel('Wavelength, nm')
ylabel('Counts')

% %% Load backgroun file
% [File, Path] = uigetfile('*.h5');   %% Insert file name
% %I_Ml_on_SiO2_0.3xname = 'PL_PAS_FIN0.7V_resolution_0.005V_0.4s_exposure_spectro750_excitation632_20µW_Olympus_002.h5_001.h5';      %% Insert file name
% name_back = strcat(Path,File);
% wvl = h5read(name_back,'/WL');
% %ev=1240/wvl;
% data_back = h5read(name,'/Data');
% background= mean(data_back);
% 
% 
% figure()
% plot(wvl,spectra-background)
% xlabel('Wavelength, nm')
% ylabel('Counts')