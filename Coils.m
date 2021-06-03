%% Use Zhenjie Yan's library to get profiles from geometries

%Note MOT is quad so anti-Helmholtz, 
%while Feshbach/bias are uniform so Helmholtz,
%and curvature ???

clear all; close all; clc;

addpath('Library/')

%% physical constants

mu=1.4e6; %1.4MHz/G 

mli=9.988346*10^-27; %kg
mer=2.7774008*10^-25; %kg

hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant

rho_c=8.96;C=0.385; %density and specific heat of copper for coils

%All units below in Gauss and cm
% K: magnetic slope/A
% G: bias field/A


%% Coil MOT (Kimball chamber)

%Geometry from Vilas (MS thesis Table 3.2), check if compatible with
%chamber geometries for mounting/optical access etc.

xMH=0;yMH=0;zMH=0; %trap location
RinMH=4.4;RoutMH=5.21;HMH=1.04; %coil geometry
nMH=[0,0,1]; %coil orientation
SepMH=3.66; %bottom-to-atom separation
DMH=SepMH+(HMH/2); %distance from atoms to center of coil

SMHa=0.1296;SMHb=0.0768; %hollow-core cross-sections, Sa-Sb = area of copper
NcoilMH=(RoutMH-RinMH)*HMH/SMHa; %number of coils
RMH=(1.7*10^(-6)*pi()*(RoutMH^2-RinMH^2)*HMH/(SMHa*SMHb)); %resistance
VMH=pi()*(RoutMH^2-RinMH^2)*HMH*SMHb/SMHa; %volume of copper

KMH=B_field_Slope_antiHelmholtz_z(0,yMH,0,xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
GMH=B_field_pair(0,yMH,0,xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
HRMH=1000/(rho_c*VMH*C); %heat capacity


%Plot z profiles
figure;

subplot(1,2,1)
z=linspace(-5,5);
GMz=z*0;
for i=1:length(z)
    temp=multi_antiHelmholtz_B(0,yMH,z(i),xMH,yMH,zMH,nMH,NcoilMH,RinMH,RoutMH,HMH,DMH);
    GMz(i)=temp(3);
end

h=plot(z*10,GMz,'Displayname','Magnetic field','linewidth',2);
hold on
%line([-35.306,-35.306],[-1,1],'color','k','Displayname','Chamber edge')
%line([35.306,35.306],[-1,1],'color','k','Displayname','Chamber edge')
line([-15,-15],[-1,1],'color','r','Displayname','MOT beam edge')
line([15,15],[-1,1],'color','r','Displayname','MOT beam edge')
hold off
title('Kimball MOT Coils');
legend('show','location','northwest')
ylabel('Magnetic field with unit current (G/A)');
xlabel('z position (mm)');xlim([-50,50])

subplot(1,2,2)
KMz = FiniteD( z,z*0,GMz,z*0,1);
h=plot(z*10,KMz,'Displayname','Magnetic field gradient','linewidth',2);
hold on
%line([-35.306,-35.306],[-1,1],'color','k','Displayname','Chamber edge')
%line([35.306,35.306],[-1,1],'color','k','Displayname','Chamber edge')
line([-15,-15],[-1,1],'color','r','Displayname','MOT beam edge')
line([15,15],[-1,1],'color','r','Displayname','MOT beam edge')
hold off
title('beam size based on Innsbruck');
legend('show','location','southwest')
ylabel('Magnetic field Gradient with unit current (G/cm/A)');ylim([-0.2,0.4])
xlabel('z position (mm)');xlim([-50,50])


%% Coil MOT_bias (Kimball chamber)

%Geometry roughly from Vilas (MS thesis Table 3.2)

xb=0;yb=0;zb=0;
Rinb=6.68;Routb=6.95;Hb=0.52;
nb=[0,0,1];
Sepb=3.66;
Db=Sepb+(Hb/2);

Sba=0.0384;Sbb=0.03;
Ncoilb=(Routb-Rinb)*Hb/Sba;
Rb=(1.7*10^(-6)*pi()*(Routb^2-Rinb^2)*Hb/(Sba*Sbb));
Vb=pi()*(Routb^2-Rinb^2)*Hb*Sbb/Sba;

Kb=B_field_Slope_antiHelmholtz_z(0,yb,0,xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
Gb=B_field_pair(0,yb,0,xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
HRb=1000/(rho_c*Vb*C);

%Plot z profiles
figure;

subplot(1,3,1)
z=linspace(-5,5);
Gbz=z*0;
for i=1:length(z)
    temp=B_field_pair(0,yb,z(i),xb,yb,zb,nb,Ncoilb,Rinb,Routb,Hb,Db);
    Gbz(i)=temp(3);
end

h=plot(z*10,Gbz,'color','b','Displayname','Magnetic field','linewidth',2);
hold on
line([-35.306,-35.306],[-1,1],'color','k','Displayname','Chamber edge')
line([35.306,35.306],[-1,1],'color','k','Displayname','Chamber edge')
hold off
legend('show','location','southwest')
ylabel('Magnetic field with unit current (G/A)');ylim([0,0.6])
xlabel('z position (mm)');xlim([-50,50])

subplot(1,3,2)
Kbz = FiniteD( z,z*0,Gbz,z*0,1);
h=plot(z*10,Kbz,'color','g','Displayname','Magnetic field gradient','linewidth',2);
hold on
line([-35.306,-35.306],[-1,1],'color','k','Displayname','Chamber edge')
line([35.306,35.306],[-1,1],'color','k','Displayname','Chamber edge')
hold off
title('Kimball Chamber Bias Coils');
legend('show','location','southwest')
ylabel('Gradient with unit current (G/cm/A)');ylim([-0.5,0.5])
xlabel('z position (mm)');xlim([-50,50])

subplot(1,3,3)
Cbz = FiniteD( z,z*0,Kbz,z*0,1);
h=plot(z*10,Cbz,'color','r','Displayname','Magnetic field curvature','linewidth',2);
hold on
line([-35.306,-35.306],[-1,1],'color','k','Displayname','Chamber edge')
line([35.306,35.306],[-1,1],'color','k','Displayname','Chamber edge')
hold off
legend('show','location','southwest')
ylabel('Curvature with unit current (G/cm/A)');ylim([-0.5,0.5])
xlabel('z position (mm)');xlim([-50,50])


%% Coil big Feshbach (Glass cell)

%LOOKING GOOD

xbF=0;ybF=24.4355;zbF=0; %DISTANCE FROM MOT CENTER TO GLASS CENTER: 24.4355 cm
%RinbF=3.38;RoutbF=5.18;HbF=0.72;
RinbF=4.18;RoutbF=5.98;HbF=0.864;
nbF=[0,0,1];
SepbF=2.04; %Chosen based on cell thickness and 5mm offset
DbF=SepbF+HbF/2;

SbFa=0.1296;SbFb=0.0768;
NcoilbF=(RoutbF-RinbF)*HbF/SbFa;

RbF=(1.7*10^(-6)*pi()*(RoutbF^2-RinbF^2)*HbF/(SbFa*SbFb));
VbF=pi()*(RoutbF^2-RinbF^2)*HbF*SbFb/SbFa;


%MAKE SURE TO CHOOSE CORRECTLY, ANTI OR PAIR
KbF=B_field_Slope_antiHelmholtz_z(0,ybF,0,xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF); 
GbF=B_field_pair(0,ybF,0,xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF);

%Plot z profiles
figure;

subplot(1,2,1)
z=linspace(-3,3);
GFz=z*0;
for i=1:length(z)
    temp=B_field_pair(0,ybF,z(i),xbF,ybF,zbF,nbF,NcoilbF,RinbF,RoutbF,HbF,DbF);
    GFz(i)=temp(3);
end

h=plot(z*10,GFz,'Displayname','Magnetic field','linewidth',2);
hold on
line([-8.955,-8.955],[-1,2*max(GFz)],'color','k','Displayname','Cell edge')
line([8.955,8.955],[-1,2*max(GFz)],'color','k','Displayname','Cell edge')
line([-30,30],[2,2],'color','r','Displayname','1000 G mark')
hold off
title('Science Cell Feshbach Coils');
legend('show','location','southwest')
ylabel('Magnetic field with unit current (G/A)');ylim([1.5,2.5])
xlabel('z position (mm)');xlim([-30,30])

subplot(1,2,2)
KFz = FiniteD( z,z*0,GFz,z*0,1);
h=plot(z*10,KFz,'Displayname','Magnetic field gradient','linewidth',2);
hold on
line([-8.955,-8.955],[4*min(KFz),4*max(KFz)],'color','k','Displayname','Cell edge')
line([8.955,8.955],[4*min(KFz),4*max(KFz)],'color','k','Displayname','Cell edge')
line([-30,30],[-0.025,-0.025],'color','r','Displayname','F3 bounds')
line([-30,30],[0.025,0.025],'color','r','Displayname','F3 bounds')
hold off
legend('show','location','southwest')
ylabel('Magnetic field Gradient with unit current (G/cm/A)');ylim([-0.5,0.5])
xlabel('z position (mm)');xlim([-30,30])


%% coil curvature (Glass cell)

%Provides z^2 DOF for compensation and confinement. 
%~100G/cm^2 gives ~58 Hz axial trapping frequency gives 2.5 micron width

xC=0;yC=24.4355;zC=0; %DISTANCE FROM MOT CENTER TO GLASS CENTER: 24.4355 cm
RinC=2.82;RoutC=5.12;HC=0.732; %Need RinC > 2.35 for microscope
nC=[0,0,1];
SepC=3.96; %Need to accomodate the 3.336cm total height feshbach coils, so SepC > 3.336
DC=SepC+HC/2;
SCa=0.1296;SCb=0.0768;
NcoilC=(RoutC-RinC)*HC/(SCa);

RC=(1.7*10^(-6)*pi()*(RoutC^2-RinC^2)*HC/(SCa*SCb));
VC=pi()*(RoutC^2-RinC^2)*HC;

KC=B_field_Slope_antiHelmholtz_z(0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC); 
GC=B_field_pair(0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC);

%Plot z profiles
figure;

subplot(1,2,1)
z=linspace(-3,3);
GCz=z*0;
for i=1:length(z)
    temp=B_field_pair(0,yC,z(i),xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC);
    GCz(i)=temp(3);
end

h=plot(z*10,GCz,'Displayname','Magnetic field','linewidth',2);
hold on
line([-8.955,-8.955],[-1,2*max(GCz)],'color','k','Displayname','Cell edge')
line([8.955,8.955],[-1,2*max(GCz)],'color','k','Displayname','Cell edge')
hold off
title('Science Cell Curvature Coils');
legend('show','location','southwest')
ylabel('Magnetic field with unit current (G/A)');ylim([-1,2*max(GCz)])
xlabel('z position (mm)');xlim([-30,30])

KCz = FiniteD( z,z*0,GCz,z*0,1);

subplot(1,2,2)
CCz = FiniteD( z,z*0,KCz,z*0,1);
h=plot(z*10,CCz,'Displayname','Magnetic field axial curvature','linewidth',2);
hold on
line([-8.955,-8.955],[4*min(CCz),4*max(CCz)],'color','k','Displayname','Cell edge')
line([8.955,8.955],[4*min(CCz),4*max(CCz)],'color','k','Displayname','Cell edge')
line([-30,30],[0.2,0.2],'color','r','Displayname','~60Hz axial freq at 500A')
hold off
legend('show','location','southwest')
ylabel('Magnetic field curvature with unit current (G/cm^2/A)');ylim([4*min(CCz),4*max(CCz)])
xlabel('z position (mm)');xlim([-30,30])


%% Axial Magnetic field curvature of Curv
AC1=second_derivative_norm_Bzz ( 0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC  );
%Radius Magnetic field gradient of Curv
AC2=second_derivative_norm_Bxx ( 0,yC,0,xC,yC,zC,nC,NcoilC,RinC,RoutC,HC,DC  );
T=AC1/norm(GC);

% trapping potential for curvature
omegaC1=sqrt(mu*AC1*200*hh/(0.5*mli*(1/100)^2))/(2*pi);
omegaC2=omegaC1/sqrt(2);


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
title('ODT axis magnetic field profiles');
legend('show','location','eastOutside')
ylabel('Magnetic field with unit current (G/A)')
xlabel('Position along ODT axis (mm)')
set(get(h,'parent'),'fontsize',12)
print('-clipboard','-dbitmap')

