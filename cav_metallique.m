clear all
close all
clc
%%
%% Build the structure
n_air = 1;
n_quartz = 1.5534; n_poly = 1.5916; 

n_cavity = 1;% refractive indexes
n_titanium = 2.636+1i*3.6520;
L_titanium = 2;
L_gold  = 58;
L_poly = 100;
L_cavity = 120;
L_air = 1e6;
L_quartz = 100;
lambda = linspace(400,1100,1000);
N = size(lambda,2);
T = zeros(1,N);
for ll = 1:N
    wvl = lambda(ll);
    [n_real , n_imag] = gold_dispersion_relation(wvl);
    n_gold = n_real + 1i*n_imag;

%% Create the structure

n = [n_air,n_quartz, n_gold, n_poly , n_cavity, n_gold,n_quartz, n_air];
L = [L_air,L_quartz,L_gold, L_poly, L_cavity, L_gold,L_quartz, L_air];

%% Parte Yannick

E_transmission = zeros(1,N);
E_reflection = zeros(1,N); 
lambda_0 = lambda(ll);

ii=size(n,2);
E{ii}=[1e-2;0];

for ii=size(n,2):-1:2

x{ii}=L(ii-1)-L(ii);
Etot{ii}=E{ii}(1)*exp(1i*2*pi*n(ii)/lambda_0*x{ii})+E{ii}(2)*exp(-1i*2*pi*n(ii)/lambda_0*x{ii});

r=(n(ii-1)-n(ii))/(n(ii-1)+n(ii));
t=2*n(ii-1)/(n(ii-1)+n(ii));

M=[[1/t r/t];[r/t 1/t]];
x0=L(ii-1);
k1=2*pi*n(ii-1)/lambda_0;
k2=2*pi*n(ii)/lambda_0;
P2=[[exp(1i*k2*x0),0];[0,exp(-1i*k2*x0)]];
P1=[[exp(-1i*k1*x0),0];[0,exp(1i*k1*x0)]];
%D=[[exp(-1i*2*pi*n(ii-1)/lambda*Lbis(ii-1)),0];[0,exp(-1i*2*pi*n(ii-1)/lambda*Lbis(ii-1))]]
%%%D=[[exp(-1i*2*pi*Lbis(ii-1)*la0/lambda),0];[0,exp(1i*2*pi*Lbis(ii-1)*la0/lambda)]]; %propagation libre


E{ii-1}=P1*M*E{ii};

end

ii=1;
x{ii}=linspace(0,L(ii),1000);
Etot{ii}=E{ii}(1)*exp(1i*2*pi*n(ii)/lambda_0*x{ii})+E{ii}(2)*exp(-1i*2*pi*n(ii)/lambda_0*x{ii});

E_transmission(ll) = E{size(n,2)}(1)/E{ii}(1);
E_reflection(ll) = E{ii}(2)/E{ii}(1);

% figure, plot(lambda, abs(E_transmission).^2)
% xlabel("longeur d'onde, nm")
% title('Transmission')

% figure, plot(lambda, abs(E_reflection).^2)
% xlabel("longeur d'onde, nm")
% title('Reflection')
T(ll) = abs(E_transmission(ll)).^2;
R(ll) = abs(E_reflection(ll)).^2;
A(ll) = 1-R(ll)-T(ll);

end
figure,plot(lambda, T);
xlabel("Longeur d'onde, nm")
title('Transmission')

figure, plot(lambda, R);
xlabel("Longeur d'onde, nm")
title('Reflection')

figure, plot(lambda, A);
xlabel("Longeur d'onde, nm")
title('Absorption')
