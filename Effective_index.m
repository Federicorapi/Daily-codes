clear all
close all
clc
%% Effective index
 radius = 300:100:2500;
 effective_index = [1.149 1.2356 1.2879 1.3202 1.3413 1.3558 1.3661 ...
     1.3738 1.3796 1.3841 1.3877 1.3906 1.3929 1.3949 1.3965 ...
    1.3979 1.3991 1.4001 1.401 1.4017 1.4024 1.403 1.4036];

figure(1), plot(radius/1000, effective_index,'-o')
xlabel('Radius HSQ, µm')
ylabel('Effective refractive index')
yline(1.41)

%% Quality factor
 Q = [104.7 159.25 266.99 425.13 768.25 1078.9 1399.1 1833.1 2193.6 2629.7 3296.3 4370.7 ...
     5305.5 6048.7 6603.0 7236.8 8437.6 10206 11945 14878 16250 16584 18125];

figure(2), plot(radius/1000, Q);
ylabel('Quality factor')
xlabel('Radius HSQ, µm')

%% Frequency
 frequency = [3.3427 3.345 3.3317 3.3184 3.3069 3.2995 3.2925 3.2871 3.2837 3.2802 3.2773 3.2750 3.2734 3.2717...
     3.2700 3.2692 3.2682 3.2672 3.2665 3.2658 3.2652 3.2646 3.2640];

figure(3), plot(radius/1000,frequency*10^4,'-o')
ylabel('Frequency, Hz')
xlabel('Radius, µm')

thickness = @(n) (910e-9)./(2.*n);
figure(4), plot(radius, thickness(effective_index));
xlabel('Radius HSQ, µm')
ylabel('HSQ thickness')

figure(5), plot(thickness(effective_index),frequency);
xlabel('Thickness HSQ')
title('f(h)')
spline_fh = spline(thickness(effective_index),frequency);

%% Correcting the HSQ thickness to reduce the frequency error
h_corrected = zeros(1,size(effective_index,2));                        %%% These two vectors will store the final values we are looking for
f_corrected = zeros(1,size(effective_index,2));

h0 = thickness(effective_index(1));
f0 = frequency(1);

spline_thickness_vs_frequency = spline(thickness(effective_index),frequency);

f1 = frequency(2);
h1 = thickness(effective_index(2));

h1_prime = h1 - 0.01*h1;
f1_prime = ppval(spline_thickness_vs_frequency,h1_prime);
c_linear_fh = polyfit([h1_prime, h1],[f1_prime, f1],1);

n = linspace(effective_index(1),effective_index(end),1000);

[minimum, index_tilde] = min(abs(polyval(c_linear_fh,thickness(n))-ppval(spline_thickness_vs_frequency,thickness(n))));
h_tilde = thickness(n(index_tilde));


h_corrected(1) = h0; 
h_corrected(2) = h_tilde;
f_corrected (1) = f0; 
f_corrected(2) = ppval(spline_thickness_vs_frequency,h_corrected(2));



%% Fitto l'indice di rifrazione dipendente da n

spline_radius_vs_effective_index = spline(radius,effective_index);

c_linear_hr = polyfit([radius(1), radius(2)],[h_corrected(1), h_corrected(2)],1);
r = linspace(radius(1), radius(end),1000);

for jj = 3:size(effective_index,2)
h_prime = polyval(c_linear_hr,radius(jj));

[minimum, n_prime_index] = min(abs(h_prime - thickness(n))); %%% n'
% Cerco r'
n_prime = n(n_prime_index);

[minimum ,index_r_prime] = min(abs(ppval(spline_radius_vs_effective_index,r)-n_prime));
r_prime = r(index_r_prime);  %%% r'


%%%% Faccio l'altra approssimazione lineare
c_linear_hr_2 = polyfit([r_prime radius(jj)],[thickness(n_prime) thickness(effective_index(jj))],1);

[minimum, r_tilde_index] =  min(abs(thickness(ppval(spline_radius_vs_effective_index,r))-polyval(c_linear_hr_2,r)));

h_tilde = thickness(ppval(spline_radius_vs_effective_index,r(r_tilde_index)));

h_corrected(jj) = h_tilde;

f_corrected(jj) = ppval(spline_thickness_vs_frequency,h_corrected(jj));
end

figure(6), plot(radius/1000,h_corrected,'o');
hold on
xlabel('radius')
ylabel('Thickness')

plot(radius/1000,thickness(effective_index))
legend('Correction','First approximation')
hold off

figure(7), plot(h_corrected,ppval(spline_fh,h_corrected),'o');
xlabel('HSQ thickness')
ylabel('Frequency')
hold on
plot(thickness(effective_index),frequency)
legend('Corrected','Experimental')
hold off
difference = h_corrected-thickness(effective_index)