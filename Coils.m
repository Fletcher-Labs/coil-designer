%% Use Zhenjie Yan's library to get profiles from geometries

%Note MOT is quad so anti-Helmholtz, 
%while Feshbach/bias are uniform so Helmholtz,
%and curvature ???

clear all; close all; clc;

addpath('Library/')

%% physical constants

mu=1.4e6; %1.4MHz/G DOES THIS NEED TO BE SPECIES SPECIFIC??

mli=9.988346*10^-27; %kg
mer=2.7774008*10^-25; %kg

hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant

rho_c=8.96;C=0.385; %density and specific heat of copper for coils

%All units below in Gauss and cm
% K: magnetic slope/A
% G: bias field/A

%BELOW IS A WORKING SKELETON THAT JUST NEEDS TO BE POPULATED WITH GEOMETRY

%% Coil MOT (Kimball chamber)

%Geometry from Vilas (MS thesis Table 3.2), check if compatible with
%chamber geometries for mounting/optical access etc.

xMH=0;yMH=0;zMH=0; %trap location
RinMH=4.4;RoutMH=5.21;HMH=1.04; %coil geometry
nMH=[0,0,1]; %coil orientation
SepMH=3.66; %bottom-to-atom separation
DMH=SepMH+(HMH/2);

SMHa=0.1296;SMHb=0.0768; %hollow-core cross-sections, Sa-Sb = area of copper
NcoilMH=(RoutMH-RinMH)*HMH/SMHa; %number of coils
RMH=(1.7*10^(-6)*pi()*(RoutMH^2-RinMH^2)*HMH/(SMHa*SMHb)); %resistance
VMH=pi()*(RoutMH^2-RinMH^2)*HMH*SMHb/SMHa; %volume of copper


KMH=B_field_Slope_antiHelmholtz_z(0,yMH,0,xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
GMH=B_field_pair(0,yMH,0,xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
HRMH=1000/(rho_c*VMH*C); %heat capacity


%% Coil MOT_bias (Kimball chamber)

%Geometry from Vilas (MS thesis Table 3.2), check if compatible with
%chamber geometries for mounting/optical access etc.

xb=0;yb=0;zb=0;
Rinb=6.58;Routb=6.85;Hb=0.52;
nb=[0,0,1];
Sepb=3.1;
Db=Sepb+(Hb/2);

Sba=0.0384;Sbb=0.03;
Ncoilb=(Routb-Rinb)*Hb/Sba;
Rb=(1.7*10^(-6)*pi()*(Routb^2-Rinb^2)*Hb/(Sba*Sbb));
Vb=pi()*(Routb^2-Rinb^2)*Hb*Sbb/Sba;

%MAKE SURE TO CHOOSE CORRECTLY, ANTI OR PAIR
Kb=B_field_Slope_antiHelmholtz_z(0,yb,0,xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
Gb=B_field_pair(0,yb,0,xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
HRb=1000/(rho_c*Vb*C);


%% Coil big Feshbach (Glass cell)

%Geometry from Fermi3, check if compatible with
%chamber geometries for mounting/optical access etc.

xbF=0;ybF=24.4355;zbF=0; %DISTANCE FROM MOT CENTER TO GLASS CENTER: 24.4355 cm
RinbF=3.38;RoutbF=5.18;HbF=0.72;
nbF=[0,0,1];
SepbF=1.6;
DbF=SepbF+HbF/2;
SbFa=0.1296;SbFb=0.0768;
NcoilbF=(RoutbF-RinbF)*HbF/SbFa;

RbF=(1.7*10^(-6)*pi()*(RoutbF^2-RinbF^2)*HbF/(SbFa*SbFb));
VbF=pi()*(RoutbF^2-RinbF^2)*HbF*SbFb/SbFa;

%MAKE SURE TO CHOOSE CORRECTLY, ANTI OR PAIR
KbF=B_field_Slope_antiHelmholtz_z(0,ybF,0,xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF); 
GbF=B_field_pair(0,ybF,0,xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF);


%% coil curvature (Glass cell)

%Geometry from Fermi3, check if compatible with
%chamber geometries for mounting/optical access etc.

xC=0;yC=24.4355;zC=0; %DISTANCE FROM MOT CENTER TO GLASS CENTER: 24.4355 cm
RinC=1.92;RoutC=4.32;HC=0.72;
nC=[0,0,1];
SepC=3.96;
DC=SepC+HC/2;
SCa=0.1296;SCb=0.0768;
NcoilC=(RoutC-RinC)*HC/(SCa);

RC=(1.7*10^(-6)*pi()*(RoutC^2-RinC^2)*HC/(SCa*SCb));
VC=pi()*(RoutC^2-RinC^2)*HC;

%MAKE SURE TO CHOOSE CORRECTLY, ANTI OR PAIR
KC=B_field_Slope_antiHelmholtz_z(0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC); 
GC=B_field_pair(0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC);

%% Axial Magnetic field curvature of Curv
AC1=second_derivative_norm_Bzz ( 0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC  );
%Radius Magnetic field gradient of Curv
AC2=second_derivative_norm_Bxx ( 0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC  );
T=AC1/norm(GC);

% trapping potential for curvature
omegaC1=sqrt(mu*AC1*200*hh/(0.5*mli*(1/100)^2))/(2*pi);
omegaC2=omegaC1/sqrt(2);

%% Big Feshbach field
z=linspace(-3,3);
GFz=z*0;
for i=1:length(z)
    temp=B_field_pair(0,ybF,z(i),xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF);
    GFz(i)=temp(3);
end

figure(1)
h=plot(z*10,GFz,'Displayname','Magnetic field','linewidth',2);
hold on
line([-8,-8],[-1,2*max(GFz)],'color','k','Displayname','Cell edge')
line([8,8],[-1,2*max(GFz)],'color','k','Displayname','Cell edge')
hold off
title('Science Cell Feshbach Coils');
legend('show','location','southwest')
ylabel('Magnetic field with unit current (G/A)');ylim([0,2*max(GFz)])
xlabel('z position (mm)')
set(get(h,'parent'),'fontsize',12)
print('-clipboard','-dbitmap')
%% Big Feshbach field gradient
KFz = FiniteD( z,z*0,GFz,z*0,1);

figure(2)
h=plot(z*10,KFz,'Displayname','Magnetic field gradient','linewidth',2);
hold on
line([-8,-8],[4*min(KFz),4*max(KFz)],'color','k','Displayname','Cell edge')
line([8,8],[4*min(KFz),4*max(KFz)],'color','k','Displayname','Cell edge')
hold off
title('Science Cell Feshbach Gradient');
legend('[4*min(KFz),4*max(KFz)]show','location','southwest')
ylabel('Magnetic field Gradient with unit current (G/cm/A)');
xlabel('z position (mm)')
ylim([4*min(KFz),4*max(KFz)])
set(get(h,'parent'),'fontsize',12)
print('-clipboard','-dbitmap')

%% Next steps

%Would be nice to have a cross-section plot showing all fields/gradients
%from Kimball to glass cell. Similar to Vilas fig 3.11

y=linspace(-10,35);
GMHy=y*0; %MOT
Gby=y*0; %Bias
GbFy=y*0; %big Fesh
Gcy=y*0; %Curv

for i=1:length(z)
    tempMH=multi_antiHelmholtz_B(0,y(i),0,xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
    tempb=B_field_pair(0,y(i),0,xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
    tempbF=B_field_pair(0,y(i),0,xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF);
    tempc=B_field_pair(0,y(i),0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC);
    
    GMHy(i)=tempMH(3);
    Gby(i)=tempb(3);
    GMbFy(i)=tempbF(3);
    Gcy(i)=tempc(3);
end

figure(1)
h=plot(y*10,GMHy,'Displayname','MOT','linewidth',2);
hold on
plot(y*10,Gby,'Displayname','Bias','linewidth',2);
plot(y*10,GbFy,'Displayname','Big Feshbach','linewidth',2);
plot(y*10,Gcy,'Displayname','Curvature','linewidth',2);

line([0,0],[-0.5,1.5],'color','k','Displayname','MOT center')
line([244.355,244.355],[-0.5,1.5],'color','m','Displayname','Science cell center')

hold off
title('Coils based on Cambridge and Fermi3 geometries');
legend('show','location','eastOutside')
ylabel('Magnetic field with unit current (G/A)')
xlabel('Position along ODT axis (mm)')
set(get(h,'parent'),'fontsize',12)
print('-clipboard','-dbitmap')

