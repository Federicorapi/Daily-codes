clear all
close all
clc
lambda = 700:100:1100;
%% Gold
k = [4.0624 4.9077 5.7227 6.4731 7.2433];
thickness = linspace(0,80);
alpha = (4*pi./lambda).*k;

figure(1)
hold on
for i = 1:size(k,2)
reflection = 1-exp(-alpha(i).*thickness);
plot(thickness, reflection)
plot(thickness, 1-reflection,'--','handleVisibility','off')

end
xlabel('Gold thicknesss, nm')
title('- -, Transmission. -- Reflection')
legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')

%% Silver
k = [4.8025 5.5698 6.3711 7.1155 7.8918];
thickness = linspace(0,80);
alpha = (4*pi./lambda).*k;

figure(2)
hold on
for i = 1:size(k,2)
reflection = 1-exp(-alpha(i).*thickness);
plot(thickness, reflection)
plot(thickness, 1-reflection,'--','handleVisibility','off')

end
xlabel('Silver thicknesss, nm')
title('- -, Transmission. -- Reflection')
legend('\lambda = 700','\lambda = 800','\lambda = 900','\lambda = 1000','\lambda = 1100')


