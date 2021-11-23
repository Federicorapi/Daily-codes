function [] = plot_spectra (x,y,name)

data = h5read(name,'/Data');
wvl = h5read(name,'/WL');
xy = h5read(name,'/Consigne');


spectra= data(:,x+(y-1)*21);


figure()
plot(wvl(226:1340),spectra(226:1340))
xlabel('Wavelength, nm')
ylabel('Counts')
titre=strcat('Spectra at pixel',' (',num2str(x),',',num2str(y),')');
subtitle('Puissance 67µW, 3s exposure time, 632.8 nm excitation')
title(titre)



