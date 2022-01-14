clear all 
close all
clc


%%
na = 1; nb = 1.52; nH = 2.6; nL = 1.4585; 
n_cavity = na;% refractive indexes
lambda = 800; %nm
LH = lambda/(4*nH);
LL = lambda/(4*nL);
L_cavity = lambda/(2*na);   % Air cavity 


Thicknesses = [1.27054000000000;1.74507000000000;1.03998000000000;1.21835000000000;...
    1.24900000000000;1.25126000000000;0.686350000000000;0.0852426000000000;1.04015000000000;...
    0.925410000000000;0.708470000000000;0.795390000000000;1.00279000000000;1.03490000000000;...
    0.996630000000000;1.00698000000000;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];



Materials = [-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2;-5;-2];

for i = 1 : size(Materials,1)
    
    if Materials(i) == -5
        Materials(i) = 2.6;
        Thicknesses(i) = Thicknesses(i)*LH;

        
    elseif Materials(i) == -2
        Materials(i) = 1.4585;
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

la = linspace(500,1500,501);
k_vector = 2*pi./la;


%% Plot Reflectivity
Gla = 100*abs(Multidiel_Federico(n,L.*n(2:end-1),k_vector)).^2;
figure(2), plot(la,Gla);
xlabel('Wavelength, nm')
ylabel('Reflectivity, A.U.')

%% Plot Refractive indexes profile
N_steps = size(L,2);
indice = 2;
start = 0;
figure(3),
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

xlabel('Bragg depth, nm')
ylabel('Refractive index')
xticks







