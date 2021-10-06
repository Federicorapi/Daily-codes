clear all
clc
close all
%% Opening of h5 file
name = 'map_ML_on_SiO2_001.h5';
data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
xy = h5read(name,'/Consigne');
%% Creation of one spectrum to see the background
sum_of_spectra = zeros(134,1);

for ii = 1:525
    sum_of_spectra = sum_of_spectra+data(:,ii);
end
figure(1),plot(1240./wvl,sum_of_spectra);
title('Sum of all the spectra')
xlabel('Energy, eV')
ylabel('A.U.')
%% Analysis of peak Frodo
background = 121600;
intensities_peak1 = ones(25,25);
r = 1;
for y = 1:21
    for x = 1:25
        intensities_peak1 (x,y) = max(data(19:28,x+(r-1)*25));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(2),
imagesc(intensities_peak1-background)
title('Intensity map of peak 3')
axis xy
xlim([1 21])
colorbar
%% Analysis of peak Pipino
background = 135100;
intensities_peak2 = ones(25,25);
r = 1;
for y = 1:21
    for x = 1:25
        intensities_peak2(x,y) = max(data(29:32,x+(r-1)*25));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(3),
imagesc(intensities_peak2-background)
title('Intensity map of peak 2')
axis xy
xlim([1 21])
colorbar
%% Analysis of peak Sauron
background = 89080;
intensities_peak3 = ones(25,25);
r = 1;
for y = 1:21
    for x = 1:25
        intensities_peak3(x,y) = max(data(60:63,x+(r-1)*25));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(4),
imagesc(intensities_peak3-background)
title('Intensity map of peak 1')
axis xy
xlim([1 21])
colorbar
%% Analysis of peak Sam
background = 80510;
intensities_peak4 = ones(25,25);
r = 1;
for y = 1:21
    for x = 1:25
        intensities_peak4(x,y) = max(data(14:19,x+(r-1)*25));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(5), title('Intensity map of Sam')
imagesc(intensities_peak4-background)
axis xy
title('Intensity map of 4')
xlim([1 21])
colorbar
%% Making the map
intensities = ones(25,25);
background = ones(25,25);
r = 1;
for y = 1:21
    for x = 1:25
        intensities (x,y) = max(data(1:80,x+(r-1)*25));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
        background (x,y) = mean(data(1:10,x+(r-1)*25));
    end
    r = r+1;
end
figure(6), imagesc(intensities)
axis xy
xlim([1 21])
xlabel('Pixels')
ylabel('Pixels')
colorbar
title('Total intensity map')

% for i = 6:16
% background = mean(data(1:25,i));
% figure(i), plot(wvl,data(:,i)-background)
% title(strcat('Pixel number',num2str(i)))
% ylim ([-100 inf])
% end

%% Altro scan
%% Opening of h5 file
name = 'ML_on_SIO2_map_exp1p5_pow225µW_grating300_centralwvl800_temp30K_001.h5';
data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
xy = h5read(name,'/Consigne');
%% Creation of one spectrum to see the background
sum_of_spectra = zeros(134,1);

for ii = 1:961
    sum_of_spectra = sum_of_spectra+data(:,ii);
end
figure(1),plot(sum_of_spectra,'LineWidth',1.2);
title('Sum of all the spectra')
title('Sum of all the spectra')
xlabel('Energy, eV')
ylabel('A.U.')
%% Analysis of peak 1
background = 610100;
intensities_peak1 = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities_peak1 (x,y) = max(data(13:18,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(2),
imagesc(intensities_peak1-background)
    title('Intensity map of peak 4')
axis xy
colorbar
%% Analysis of peak 2
background = 643700;
intensities_peak2 = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities_peak2(x,y) = max(data(21:26,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(3),
imagesc(intensities_peak2-background)
title('Intensity map of peak 3')
axis xy
colorbar
%% Analysis of peak 3
background = 635100;
intensities_peak3 = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities_peak3(x,y) = max(data(28:33,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(4),
imagesc(intensities_peak3-background)
title('Intensity map of peak 2')
axis xy
colorbar
%% Analysis of peak 4
background = 625900;
intensities_peak4 = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities_peak4(x,y) = max(data(79:87,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
    end
    r = r+1;
end
figure(5),
imagesc(intensities_peak4-background)
title('Intensity map of peak 1')
axis xy
colorbar
%% Making the map of the substrate
intensities = ones(31,31);
background = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities (x,y) = max(data(1:5,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
        background (x,y) = mean(data(116:130,x+(r-1)*31));
    end
    r = r+1;
end
figure(6), imagesc(intensities-background)
axis xy
xlim([1 31])
xlabel('Pixels')
ylabel('Pixels')
title('Map of the substrate')
colorbar

%% Making the total map
intensities = ones(31,31);
background = ones(31,31);
r = 1;
for y = 1:31
    for x = 1:31
        intensities (x,y) = max(data(10:134,x+(r-1)*31));
        %intensities(x,y) = trapz(data(1:50,x+(r-1)*25));
        background (x,y) = mean(data(116:130,x+(r-1)*31));
    end
    r = r+1;
end
figure(7), imagesc(intensities-background)
axis xy
xlim([1 31])
xlabel('Pixels')
ylabel('Pixels')
colorbar

title('Total intensity map')
% for i = 6:16
% background = mean(data(1:25,i));
% figure(i), plot(wvl,data(:,i)-background)
% title(strcat('Pixel number',num2str(i)))
% ylim ([-100 inf])
% end






