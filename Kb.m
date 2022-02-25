function[kwb,Kw]=Kb(Hs,Va,S,T);
% Input: 
% Hs - significant waveheight (m)
% Va - Volume of air entrainment rate by breaking (m/s)
% T  - seawater temperature (deg C)
% S  - seawater salinity (psu)
% Output:
% Kb - bubble-mediated gas transfer coefficient (m/s) 
% Kw - total transfer coefficient including the molecular-turbulent contribution (m/s)

% Kb is computed following Deike and Melville 2018 with optimization method developed
% by Shin Y, Deike L, and Romero L (manuscript in preparation)
%
% The optimization is based on a similarity relation for the the effective volume fraction (delta) 
% contributing to the gas exchange as a function of the dimensionless wave height (Hs/L). 
%
% Youngmi Shin - U. Connecticut
% L Romero - U. Connecticut


[Sc,nu,alpha]=SW_ND3(T,S);
D=nu./Sc;
L=(D.^2./9.8).^(1/3);

x1=Hs./L; % non-dimensional wave height 
x=log(x1);

%        p1 =  -0.0001479;%  (-0.000149, -0.0001469)
%        p2 =    0.007831;%  (0.007773, 0.007888)
%        p3 =     -0.1495;%  (-0.1507, -0.1483)
%        p4 =       1.239;%  (1.227, 1.25)
%        p5 =      -3.786;%  (-3.824, -3.747)
%        
%  Fmm = p1.*x.^4 + p2.*x.^3 + p3.*x.^2 + p4.*x + p5;

p11 =  -0.0001505327358;%  (-0.000149, -0.0001469)
   p22 =    0.00792478703;%  (0.007773, 0.007888)
       p33 =     -0.1503422767;%  (-0.1507, -0.1483)
       p44 =       1.2361506;%  (1.227, 1.25)
       p55 =      -3.740647923;%  (-3.824, -3.747)
       
 Fmm= p11.*x.^4 + p22.*x.^3 + p33.*x.^2 + p44.*x + p55;

%  f=find(x<9.5); % when log(Hs*) is less than 9.5 => Fmm should be 0.
%  Fmm(f)=0;
 Fmm=max(Fmm,0);

% a =     0.02794;%  (0.02781, 0.02808)
% b =      0.1497;%  (0.1494, 0.1499)
% c =     -0.1652;%  (-0.1656, -0.1648)
% 
%        
%  Fmm= a.*x.^b+c;% power law
% 
% [f1,f2]=find(Fmm<0);
% mm=min(abs(Fmm(:)));
% Fmm(f1,f2)=abs(mm);% remove the minus value

F_model=Fmm;

un=360000;
     Kbout1=Va.*F_model.*coeff./(2.*pi.*alpha);%m/s
     kwb=Kbout1.*sqrt(Sc/660)*un;
% total gas transfer coefficient including molecular-tubulent contribution (Deike and Melville 2018; Fairall et al. 2011)
      Anb=1.55*10^-4;
      Kw660=Anb.*Ust+Kbout1.*sqrt(Sc./660); 
      Kw=(Kw660).*un;
end
