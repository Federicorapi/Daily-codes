clear all
clc
close all
%% Opening of h5 file
name = 'LASER_ML_on_SiO2_20�W_0.004exposure_spectro632_excitation_632nm_AVEC_WINDOW_FOCALISE_001.h5';      %% Insert file name
data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
ev=1240/wvl;
xy = h5read(name,'/Consigne');
pixel=length(wvl);
scansizex=81;       %% Insert number of steps along x
scansizey=51;       %% Insert number of steps along y

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
        intensities (x,y) = max(data(1:pixel,x+(r-1)*scansizex));
        intensitiestrapz(x,y) = trapz(data(1:pixel,x+(r-1)*scansizex));
        %background (x,y) = mean(data(1:17,x+(r-1)*S));
    end
    r = r+1;
end
figure(4)
imagesc(intensities)
%% Muon killing

treshold=3000;
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

figure(5)
imagesc(intensities)
title('Map with muons killed')


%% Normalisation
%
% normfactor=21e-6*0.4;   %% change this according to excitation power and exposure time
% intensities=intensities./normfactor;
% intensities=intensities./max(max(intensities));
% data=data./normfactor;
% data=data./max(max(data));

%% Plot of a slice

tranche=zeros(1,scansizex);
yline =5;
xline =0;

for jj=1:scansizex
    tranche(jj) = intensities(jj,yline);
end

% for jj=1:scansizey
%     tranche(jj) = intensities(xline,jj);
% end


%m=mean(tranche(40:81));
m=min(tranche);
%tranche=(tranche-m);
%M=mean(tranche(1:28));
M=max(tranche);
%tranche=tranche./M;


% tranche(40:81)=m;
% tranche=(tranche-m);
% tranche(1:28)=M;
% tranche=tranche./M;
lignex=2*2.64*3*linspace(1,scansizex,scansizex);
colonney=0.7*2*2.64*3*linspace(1,scansizey,scansizey);


%tranchecut=tranche(1:20);
%lignexcut=lignex(1:20);

figure(6)
% plot(colonney,tranche)
 plot(lignex,tranche)
%plot(lignexcut,tranchecut)

%c*(1/2)*(1+erf((x-a)/(sqrt(2)*b)))

% load('fitnowindow.mat')
% load('fitwithwindow.mat')
% figure()
% plot(fitnowindow)
% title('Fit no window, on ML')
% xticks('auto')
% yticks('auto')
% xlabel('Microns')
% ylabel('Microns')
% figure()
% plot(fitwithwindow,'b')
% title('Fit with window, on ML')
% xticks('auto')
% yticks('auto')
% xlabel('Microns')
% ylabel('Microns')

%% Rescaling of the axis

xy(1,:)=2*2.64*3*xy(1,:);
xy(2,:)=0.7*2*2.64*3*xy(2,:);

% Plot of the max intensity map
figure(7),
imagesc(transpose(intensities))
%imagesc(xy(1,:),xy(2,:),transpose(intensities))
axis xy
axis equal
xticks('auto')
yticks('auto')
xlabel('pixels')
ylabel('pixels')

colorbar
title('Laser scan of ML on SiO2. 20 \muW, 0.004 seconds. Olympus objective')
ylim([1 scansizey-1])

%% Region of interest on the map
%% Select the pixel interval to create the ROI
limitx1=41;
limitx2=46;
limity1=3;
limity2=8;
spectraROI= 0;

    for x = limitx1:limitx2
        for y = limity1:limity2
            spectraROI=spectraROI+intensitiestrapz(x,y);
        end
    end
