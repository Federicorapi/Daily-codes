function [Gamma,Z]= Multidiel_Federico(n,L,k_vector,theta,pol)



% na = 1; nb = 1.52; nH = 2.076; nL = 1.468;
% LH = 0.25;
% LL = 0.25;
% N = 8;
% n = [na, nH, repmat([nL,nH], 1, N), nb]; 
% L = [LH, repmat([LL,LH], 1, N)];




if nargin == 0, print("Not enough input arguments! Session closed"); return; end   % Manage the input variables
if nargin <= 4, pol = "te"; end
if nargin == 3, theta = 0; end

if size(n,2) == 1, n = n'; end                                                     % Adjust n shape if needed

M = size(n,2)-2;

theta = theta*pi/180;


if pol == "te"
    Nsin2 = (n(1)*sin(theta))^2;                                                 % Nsin2 is a constant value (from Snell law n1sin1 = n2sin2). The square is just a conveniecy
    c = sqrt(1-Nsin2./n.^2);                                                  % cos(theta) 
    nT = n.*c;
    r = n2r(nT);                                                                   % Fresnel
else
    Nsin2 = (n(1)^2 * sin(theta))^2 / (n(1)^2 * cos(theta)^2 + n(1)^2 * sin(theta)^2);
    c = sqrt(1 - Nsin2 ./ n.^2);
    nTinv = c ./ n;                                                           % nTinv(i) = 1/nT(i) to avoid NaNs at theta = 90°
    r = -n2r(nTinv);                                                               % minus sign because n2r(n) = -n2r(1./n) 
end


L = L .*c(2:M+1);                                                              % Polarisation dependent optical lengths


Gamma = r(M+1) * ones(1,length(k_vector));

for i = M:-1:1
delta = L(i).*k_vector;
z = exp(2*1j*delta);
Gamma = (r(i)+Gamma.*z)./(1+r(i)*Gamma.*z);
end

Z = (1+Gamma)./(1-Gamma);
end





