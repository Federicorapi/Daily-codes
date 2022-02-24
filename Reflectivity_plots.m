clear all
close all
clc
lambda = 700:100:1100;
%% Gold
%% Reflection
n = ([0.131+1i*4.0624 0.15352+1i*4.9077 0.17435+1i*5.7227 0.22769+1i*6.4731 0.27750+1i*7.2433]);
thickness = linspace(0,80);
alpha = (2*pi./lambda);
nair = 1;
r_tot = zeros(size(n,2),100);
for ii = 1:size(lambda,2)
    for j = 1:100
        r11 = (nair-n(ii))/(nair+n(ii));
        r22 = (n(ii)-nair)/(nair+n(ii));
        t12 = 2*nair/(n(ii)+nair);
        t21 = 2*n(ii)/(n(ii)+nair);
        loop=r22^2*exp(2*1i*alpha(ii)*n(ii)*thickness(j));
        r_tot(ii,j) = r11+(t12*t21*r22*exp(2*1i*alpha(ii)*n(ii)*thickness(j)))*1/(1-loop);
        j
    end
   ii
end

R = abs(r_tot).^2;
figure(1),
for i = 1:5
    plot(thickness,R(i,:))
    hold on
end
title('Reflection')
xlabel('Thickness, nm')
legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')

% Transmission
t_tot = zeros(size(n,2),100);
for ii = 1:size(lambda,2)
    for j = 1:100
        r11 = (nair-n(ii))/(nair+n(ii));
        r22 = (-nair+n(ii))/(nair+n(ii));
        t12 = 2*nair/(n(ii)+nair);
        t21 = 2*n(ii)/(n(ii)+nair);
        loop=r22^2*exp(2*1i*alpha(ii)*n(ii)*thickness(j));
        t_tot(ii,j) = t12*t21*exp(1i*alpha(ii)*n(ii)*thickness(j))*(1/(1-loop));
        j
    end
   ii
end
T = abs(t_tot).^2;
A =1-( abs(t_tot).^2+abs(r_tot).^2);
figure(2),
for i = 1:5
    plot(thickness,T(i,:))
    hold on
end
title('Transmission')
xlabel('Thickness, nm')

legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')


figure(3),
for i = 1:5
    plot(thickness,A(i,:))
    hold on
end
title('Absorbance')
xlabel('Thickness, nm')

legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')


%% Silver
n = [0.041000+1i*4.8025 0.036759+1i*5.5698 0.040000+1i*6.3711 0.040000+1i*7.1155  0.044688+1i*7.8918];
%% Reflection
thickness = linspace(0,80);
alpha = (2*pi./lambda);
nair = 1;
r_tot = zeros(size(n,2),100);
for ii = 1:size(lambda,2)
    for j = 1:100
        r11 = (nair-n(ii))/(nair+n(ii));
        r22 = (n(ii)-nair)/(nair+n(ii));
        t12 = 2*nair/(n(ii)+nair);
        t21 = 2*n(ii)/(n(ii)+nair);
        loop=r22^2*exp(2*1i*alpha(ii)*n(ii)*thickness(j));
        r_tot(ii,j) = r11+(t12*t21*r22*exp(2*1i*alpha(ii)*n(ii)*thickness(j)))*1/(1-loop);
        j
    end
   ii
end

R = abs(r_tot).^2;
figure(4),
for i = 1:5
    plot(thickness,R(i,:))
    hold on
end
title('Reflection')
xlabel('Thickness, nm')

legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')

% Transmission
t_tot = zeros(size(n,2),100);
for ii = 1:size(lambda,2)
    for j = 1:100
        r11 = (nair-n(ii))/(nair+n(ii));
        r22 = (-nair+n(ii))/(nair+n(ii));
        t12 = 2*nair/(n(ii)+nair);
        t21 = 2*n(ii)/(n(ii)+nair);
        loop=r22^2*exp(2*1i*alpha(ii)*n(ii)*thickness(j));
        t_tot(ii,j) = t12*t21*exp(1i*alpha(ii)*n(ii)*thickness(j))*(1/(1-loop));
        j
    end
   ii
end
T = abs(t_tot).^2;
A =1-( abs(t_tot).^2+abs(r_tot).^2);
figure(5),
for i = 1:5
    plot(thickness,T(i,:))
    hold on
end
title('Transmission')
xlabel('Thickness, nm')

legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')


figure(6),
for i = 1:5
    plot(thickness,A(i,:))
    hold on
end
title('Absorbance')
xlabel('Thickness, nm')

legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')





