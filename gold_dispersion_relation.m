function [v_real, v_imag] = gold_dispersion_relation (wvl)
 %% Copy - past donnés discretes
lambda = [0.1879
0.1916
0.1953
0.1993
0.2033
0.2073
0.2119
0.2164
0.2214
0.2262
0.2313
0.2371
0.2426
0.249
0.2551
0.2616
0.2689
0.2761
0.2844
0.2924
0.3009
0.3107
0.3204
0.3315
0.3425
0.3542
0.3679
0.3815
0.3974
0.4133
0.4305
0.4509
0.4714
0.4959
0.5209
0.5486
0.5821
0.6168
0.6595
0.7045
0.756
0.8211
0.892
0.984
1.088
1.216
1.393
1.61
1.937];

n = [1.28
1.32
1.34
1.33
1.33
1.3
1.3
1.3
1.3
1.31
1.3
1.32
1.32
1.33
1.33
1.35
1.38
1.43
1.47
1.49
1.53
1.53
1.54
1.48
1.48
1.5
1.48
1.46
1.47
1.46
1.45
1.38
1.31
1.04
0.62
0.43
0.29
0.21
0.14
0.13
0.14
0.16
0.17
0.22
0.27
0.35
0.43
0.56
0.92];

k = [1.188
1.203
1.226
1.251
1.277
1.304
1.35
1.387
1.427
1.46
1.497
1.536
1.577
1.631
1.688
1.749
1.803
1.847
1.869
1.878
1.889
1.893
1.898
1.883
1.871
1.866
1.895
1.933
1.952
1.958
1.948
1.914
1.849
1.833
2.081
2.455
2.863
3.272
3.697
4.103
4.542
5.083
5.663
6.35
7.15
8.145
9.519
11.21
13.78
];
k = 1*i*k;
n = n+k;
n_gold = n.';
lambda = lambda.';
lambda = lambda*1e3;

%% Plot discrete values
% figure(777),
% plot(lambda, real(n_gold))
% hold on 
% plot(lambda, imag(n_gold))
% legend('n','k');
% xlabel('Wavelength,nm')
% title('Gold refractive index')


%% Fit
%lambda_fit = linspace(lambda(1),lambda(end),1000);
s_real = spline(lambda,real(n_gold));
v_real = ppval(s_real,wvl);
% figure, plot(wvl,v_real)
% hold on 
% plot(lambda,real(n_gold))
% legend('Fit','Data')
% title('Piecewise fit, real part')


%lambda_fit = linspace(lambda(1),lambda(end),1000);
s_imag = spline(lambda,imag(n_gold));
v_imag = ppval(s_imag,wvl);
% figure, plot(wvl,v_imag)
% hold on 
% plot(lambda,imag(n_gold))
% legend('Fit','Data')
% title('Piecewise fit, imaginary part')


