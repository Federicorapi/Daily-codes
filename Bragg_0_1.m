clear all 
close all
clc
%% Build the structure
na = 1; nb = 1.52; nH = 2.34; nL = 1.477; 
n_cavity = na;% refractive indexes
lambda_thickness = 1386; %nm
LH = lambda_thickness/(4*nH);
LL = lambda_thickness/(4*nL);
L_cavity = lambda_thickness/(2*na);   % Air cavity 



Thicknesses = [];

Materials = [];

for i = 1 : size(Materials,1)
    
    if Materials(i) == -5
        Materials(i) = nH;
        Thicknesses(i) = Thicknesses(i)*LH;

        
    elseif Materials(i) == -2
        Materials(i) = nL;
        Thicknesses(i) = Thicknesses(i)*LL;
    end
end
Thicknesses = Thicknesses';
Thicknesses = flip(Thicknesses,2);
Materials = Materials';
Materials = flip(Materials, 2);


N = size(Materials,2);
n = [na, Materials , n_cavity, flip(Materials), na]; 
L = [Thicknesses, L_cavity, flip(Thicknesses)]; 

% n = [na, Materials,nb]; 
% L = [Thicknesses]; 

la = linspace(1000, 1650,5001);
k_vector = 2*pi./la;
%% Plot Refractive indexes profile
N_steps = size(L,2);
indice = 2;
start = 0;
figure(2),
plot(zeros(1,100),linspace(n(1),n(2),100),'b')
hold on
for i = 1:N_steps
    Bragg_depth = sum(L(1:i),2);       %%%% Horizontal lines
    plot(linspace(start,Bragg_depth,100),ones(1,100)*n(indice),'b')
    start = Bragg_depth;
    indice = indice +1;
    hold on
end


indice = 2;
for i = 1:N_steps 
    Bragg_depth = sum(L(1:i),2);
    plot(ones(1,100)*Bragg_depth,linspace(n(indice),n(indice+1),100),'b')
    indice = indice +1;
end
hold off

xlabel('Bragg depth, nm')
ylabel('Refractive index')
%% Plot Reflectivity over one angle
Gla = abs(Multidiel_Federico(n,L.*n(2:end-1),k_vector,0,'te')).^2;
figure, plot(la,1-Gla)
xlabel('Wavelength, nm')
ylabel('Trabsmission, A.U.')
title('TE mode')

Gla = abs(Multidiel_Federico(n,L.*n(2:end-1),k_vector,0,'tm')).^2;
figure,plot(la,1-Gla)
xlabel('Wavelength, nm')
ylabel('Transmission, A.U.')
title('TM mode')
%% Plot Reflectivity over all angles
minimal_angle = 0;
maximum_angle = 90;
angles = minimal_angle:1:maximum_angle;
plots_TE = zeros(size(angles,2),size(la,2));
row = 1;
for theta = minimal_angle:maximum_angle
Gla = abs(Multidiel_Federico(n,L.*n(2:end-1),k_vector,theta,'te')).^2;
plots_TE(row,:) = Gla;
row = row +1;
end

[X,Y] = meshgrid(angles, la);
figure, surf(X,Y, abs(1-plots_TE)','edgecolor','none');
view (2);
colorbar

xlabel('Degrees')
ylabel('Wavelength, nm')
zlabel('Reflectivity, A.U.')
title('TE mode')


% new_colors = zeros(size(angles,2),3);
% for kk =1 : size(angles,2)
%     newcolors(kk,:) = [kk/size(angles,2) 0 1-kk/size(angles,2)];
% end
% colororder(newcolors) 



plots_TM = zeros(size(angles,2),size(la,2));
row = 1;

for theta = minimal_angle:maximum_angle
Gla = abs(Multidiel_Federico(n,L.*n(2:end-1),k_vector,theta,'tm')).^2;
plots_TM(row,:) = Gla;

row = row+1;
end
figure, surf(X,Y, abs(1-plots_TM)','edgecolor','none');
view (2);
colorbar


xlabel('Degrees')
ylabel('Wavelength, nm')
zlabel('Reflectivity, A.U.')
title('TM mode')
% for kk =1 : size(angles,2)
%     newcolors(kk,:) = [kk/size(angles,2) 0 1-kk/size(angles,2)];
% end
% colororder(newcolors)
%% Follow up of the maximum
data = 1-plots_TE(1,:);

[pks,locs_index] = findpeaks(data);
[pks,locs_wvl] = findpeaks(data,la);

disp('Peaks are at:')
disp(locs_wvl)
wvl = input('Which one would you like to follow? \n');
k = find(la == locs_wvl(wvl));
peaks = zeros(1,size(angles,2));
peaks(1) = locs_wvl(wvl);
for theta = 2:maximum_angle
    data = 1-plots_TE(theta,:);
    [pks,locs_wvl] = findpeaks(data,la);
    wvl_diff = abs(locs_wvl-peaks(theta-1));
    [M,I] = min(wvl_diff);
    peaks(theta) = locs_wvl(I);
end
figure, plot(angles(1:end-1),2*pi.*peaks(1:end-1)/3e8)
xlabel('Angles, deg')
ylabel('Frequency,\omega=2pi\lambda/c')
hold on




