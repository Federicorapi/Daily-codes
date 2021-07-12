clear all
close all
clc
%% Antenna L3C1 2600 uC/cm^2

x = 0;
for i = 1:19
str = strcat('2021-06-10/L3C1/0_',num2str(x),'.sif');
A(i) = sifreadnk(str);
x = x + 20;
back = mean(A(i).imageData(1:536));
Z_45_B(i) = trapz(A(i).imageData(704:988)-back);
end

x = 20;
for i = 2:18
str = strcat('2021-06-10/L3C1/0_',num2str(x),'_bis.sif');
A(i) = sifreadnk(str);
x = x + 20;
back = mean(A(i).imageData(1:536));
Z_45_B_bis(i) = trapz(A(i).imageData(704:988)-back);
end
Z_45_B_bis(1) = Z_45_B(1);
Z_45_B_bis(19) = Z_45_B(19);

angles = [0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360];
angles_rad = angles*pi/180;


angles_rad(20:38) = angles_rad;
Z = zeros(1,38); 
Z(1:19) = Z_45_B;
Z(20:38) = Z_45_B_bis;
Z = Z/max(Z);

ff = @(para, x) para(1)*(cos(x-para(2))).^2 + para(3);
para_0 = [0.1908 pi/6 0.7464];
coeff = nlinfit(angles_rad,Z,ff,para_0);

figure, polarplot(angles_rad,Z,'o','MarkerFaceColor','b','LineWidth',1)
hold on
polarplot(0:0.1:2*pi,ff(coeff,0:0.1:2*pi),'r','LineWidth',1.2);
title('Nano-bridge, 45° excitation')

%% Antenna L3C1 2400 uC/cm^2
% 45
x = 0;
for i = 1:19
str = strcat('2021-06-10/L3C1_2400/0_',num2str(x),'.sif');
A(i) = sifreadnk(str);
x = x + 20;
back = mean(A(i).imageData(1:406));
Z_45_G(i) = trapz(A(i).imageData(722:1081)-back);
end
%%% Ritorno
x = 20;
for i = 2:18
str = strcat('2021-06-10/L3C1_2400/0_',num2str(x),'_bis.sif');
A(i) = sifreadnk(str);
x = x + 20;
back = mean(A(i).imageData(1:536));
Z_45_G_bis(i) = trapz(A(i).imageData(704:988)-back);
end
Z_45_G_bis(1) = Z_45_G(1);
Z_45_G_bis(19) = Z_45_G(19);
angles = [0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360];


angles_rad = angles*pi/180;
angles_rad(20:38) = angles_rad;
Z = zeros(1,38); 
Z(1:19) = Z_45_G;
Z(20:38) = Z_45_G_bis;
Z = Z/max(Z);
ff = @(para, x) para(1)*(cos(x-para(2))).^2 + para(3);
para_0 = [0.3491 pi/6 0.5809];
coeff = nlinfit(angles_rad,Z,ff,para_0);

figure, polarplot(angles_rad,Z,'o','MarkerFaceColor','b','LineWidth',1)
hold on
polarplot(0:0.1:2*pi,ff(coeff,0:0.1:2*pi),'r','LineWidth',1.2);
title('Nano-gap, 45° excitation')









