clear all
clc
close all
%% Opening of h5 file
name = 'scan_gqd_67µW_300t_3sec_001.h5';
data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
ev=1240/wvl;
xy = h5read(name,'/Consigne');

%% Creation of one spectrum to see the background
sum_of_spectra = zeros(1340,1);

for ii = 1:441
    sum_of_spectra = sum_of_spectra+data(:,ii);
end
figure(1),plot(wvl,sum_of_spectra);
title('Sum of all the spectra')
xlabel('Wavelength, nm')
ylabel('A.U.')
subtitle('Puissance 67µW, 3s exposure time, 632.8 nm excitation')

%% Making the map
intensities = ones(21,21);
background = ones(21,21);
r = 1;
for y = 1:21
    for x = 1:21
        intensities (x,y) = max(data(226:1340,x+(r-1)*21));
        intensitiestrapz(x,y) = trapz(data(226:1340,x+(r-1)*21));
        background (x,y) = mean(data(1:170,x+(r-1)*21));
    end
    r = r+1;
end
figure(6), imagesc(intensities-background)
axis xy
xlim([1 21])
xlabel('Pixels')
ylabel('Pixels')
colorbar
subtitle('Puissance 67µW, 3s exposure time, 632.8 nm excitation')
title('Max intensity map')

figure(7), imagesc(intensitiestrapz-background)
axis xy
xlim([1 21])
xlabel('Pixels')
ylabel('Pixels')
colorbar
subtitle('Puissance 67µW, 3s exposure time, 632.8 nm excitation')
title('Integrated intensity map (trapz)')



