function[kwb,Kw]=Kb(Hs,Va,T);
% Input: 
% Hs - significant wave height (m) or effective breaking height (m)
% Va - Volume of air entrainment rate by breaking (m/s)
% T  - seawater temperature (deg C)
% Output:
% Kb - bubble-mediated gas transfer coefficient (m/s) 
% Kw - total transfer coefficient including the molecular-turbulent contribution (m/s)

% Kb is computed following Deike and Melville 2018 with optimization method developed
% by Shin Y, Deike L, and Romero L (submitted to GRL)
%
% The optimization is based on a similarity relation for the the effective volume fraction (delta) 
% contributing to the gas exchange as a function of the dimensionless wave height (Hs/L)
% 
% L= (D.^2./g).^(1/3)*(Sc/alpha).^(1/3); g is gravity, D is the molecular
% diffusivity, Sc is the Schmidt number, and alpha the Ostwald dimensionless
% solubility.
%
% S= 35 psu.
%
% Youngmi Shin - U. Connecticut
% L Romero - U. Connecticut

S=35;
g=9.81;

[sol,alpha,Sc] = solco2(T);% COARE algorithm (Sc - Wanninkhov 1992;  for S=35
nu = SW_Kviscosity2(T,S);%
%solco2.m and SW_Kviscosity2.m can be downloaded from NOAA (ftp://ftp1.esrl.noaa.gov/BLO/Air-Sea/bulkalg/cor3_6/)

D=nu./Sc;% modeluclar diffusivity
L=(D.^2./g).^(1/3)*(Sc/alpha).^(1/3); % Length scale

x1=Hs./L; % non-dimensional wave height 
x=log(x1);

psc = [  3.3296e-05  -1.9139e-3   4.3052e-2  -5.2715e-1   4.3196 -20.0518];

Delta=exp(psc(1)*x.^5+psc(2)*x.^4+psc(3)*x.^3+psc(4)*x.^2+psc(5)*x.^1+psc(6));

kwb=Va.*Delta.*1./(alpha);% bubble-mediated gas transfer (cm/hr)

% total gas transfer coefficient including molecular-tubulent contribution (Deike and Melville 2018; Fairall et al. 2011)
kmt=1.55*10^-4.*sqrt(660./Sc); %  molecular-tubulent contribution
Kw=kmt+kwb; % total
end
