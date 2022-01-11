clear all
clc
close all
%% Opening of h5 file
[File, Path] = uigetfile('*.h5');   %% Insert file name
%I_Ml_on_SiO2_0.3xname = 'PL_PAS_FIN0.7V_resolution_0.005V_0.4s_exposure_spectro750_excitation632_20µW_Olympus_002.h5_001.h5';      %% Insert file name
name = strcat(Path,File);
data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
ev=1240/wvl;
xy = h5read(name,'/Consigne');
pixel=length(wvl);
step_volt = xy(1,2)-xy(1,1);
scansizex=int64((xy(1,end)-xy(1,1))/step_volt+1);       %% Insert number of steps along x
scansizey=int64((xy(2,end)-xy(2,1))/step_volt+1);      %% Insert number of steps along y

%% Creation of one spectrum to see the background
sum_of_spectra = zeros(pixel,1);

for ii = 1:length(xy(1,:))
    sum_of_spectra = sum_of_spectra+data(:,ii);
end
figure(1),plot(wvl,sum_of_spectra);
title('Sum of all the spectra')
xlabel('Wavelength, nm') 
ylabel('A.U.')


%% Making the map
intensities = ones(scansizex,scansizey);
background = zeros(scansizex,scansizey);
r = 1;
for y = 1:scansizey
    for x = 1:scansizex
        intensities(x,y) = max(data(1:pixel,x+(r-1)*scansizex));
        intensitiestrapz(x,y) = trapz(data(1:pixel,x+(r-1)*scansizex));
        %background (x,y) = mean(data(1:17,x+(r-1)*S));
    end
    r = r+1;
end
% figure(4)
% imagesc(intensities)
% colorbar
%% Muon killing

treshold=1300;
for ss=1:scansizey
    for kk=1:scansizex
        if max(intensities(kk,ss))>treshold
            if kk ~= 1 && kk~= scansizex && ss ~= 1 && ss ~= scansizey
            intensities(kk,ss)=mean([intensities(kk-1,ss-1),intensities(kk,ss-1),intensities(kk+1,ss-1),intensities(kk-1,ss),intensities(kk+1,ss), ...
                intensities(kk-1,ss+1),intensities(kk,ss+1),intensities(kk+1,ss+1)]);
            elseif ss ~= 1 && ss ~= scansizey && (kk == 1 || kk == scansizex)
              intensities(kk,ss) = mean([intensities(kk,ss-1),intensities(kk,ss+1)]);
            elseif kk ~= 1 && kk ~= scansizex && (ss == 1 || ss == scansizey)
              intensities(kk,ss) = mean([intensities(kk-1,ss),intensities(kk+1,ss)]);
            end
        end
    end
end
% intensities(1,end) = intensities(2,end);

% figure(5)
% imagesc(intensities)
% title('Map with muons killed')
% colorbar


%% Normalisation
%
% normfactor=21e-6*0.4;   %% change this according to excitation power and exposure time
% intensities=intensities./normfactor;
% intensities=intensities./max(max(intensities));
% data=data./normfactor;
% data=data./max(max(data));

% %% Plot of a slice
% 
% tranche=zeros(1,scansizex);
% yline =5;
% xline =0;
% 
% for jj=1:scansizex
%     tranche(jj) = intensities(jj,yline);
% end
% 
% % for jj=1:scansizey
% %     tranche(jj) = intensities(xline,jj);
% % end
% 
% 
% %m=mean(tranche(40:81));
% m=min(tranche);
% %tranche=(tranche-m);
% %M=mean(tranche(1:28));
% M=max(tranche);
% %tranche=tranche./M;
% 
% 
% % tranche(40:81)=m;
% % tranche=(tranche-m);
% % tranche(1:28)=M;
% % tranche=tranche./M;
% lignex=2*2.64*3*linspace(1,scansizex,scansizex);
% colonney=0.7*2*2.64*3*linspace(1,scansizey,scansizey);
% 
% 
% %tranchecut=tranche(1:20);
% %lignexcut=lignex(1:20);
% 
% figure(6)
% % plot(colonney,tranche)
%  plot(lignex,tranche)
% %plot(lignexcut,tranchecut)
% 
% %c*(1/2)*(1+erf((x-a)/(sqrt(2)*b)))
% 
% % load('fitnowindow.mat')
% % load('fitwithwindow.mat')
% % figure()
% % plot(fitnowindow)
% % title('Fit no window, on ML')
% % xticks('auto')
% % yticks('auto')
% % xlabel('Microns')
% % ylabel('Microns')
% % figure()
% % plot(fitwithwindow,'b')
% % title('Fit with window, on ML')
% % xticks('auto')
% % yticks('auto')
% % xlabel('Microns')
% % ylabel('Microns')
% 
%% Rescaling of the axis
objective = 1; %Put 1 for Mitutotyo, 2 for Olympus
xy(1,:)=objective*2*2.64*3*xy(1,:); %from volt to micron  %MULTIPLY BOTH BY 2 TO USE THE OLYMPUS
xy(2,:)=objective*0.7*2*2.64*3*xy(2,:);   %from volt to micron

% Plot of the max intensity map
figure(7),
% imagesc(transpose(intensities))
intensities_transpose = transpose(intensities);
imagesc(xy(1,:),xy(2,:),intensities_transpose)
axis xy
axis equal
xticks('auto')
yticks('auto')
xlabel('Micron')
ylabel('Micron')

colorbar
titolo = strcat('Scan 0.4 seconds, 20\mu W excitation,',num2str(step_volt),'V steps');
title(titolo)
ylim([xy(2,1) xy(2,end)])
xlim([xy(1,1) xy(1,end)])
figure(8), imagesc(intensities_transpose)
title('Image with pixel reference')
%% Evaluation of resolution

tranche = 38; % Inserire qui il numero di riga su cui tagliare
xdirection = [1:30];
cut = intensities_transpose(tranche,xdirection);
cut = cut-min(cut);
cut = cut./max(cut);
figure(9), plot(xdirection,cut)
