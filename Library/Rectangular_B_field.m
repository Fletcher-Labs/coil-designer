function [ Bout ] = Rectangular_B_field( x,y,z,x0,y0,z0,n0,I,W,L,T )
% Calculates the magnetic field due to a rectangular pair of coils
%   XYZ - point of calculation
%   x0y0z0 - coil pair center
%   n0 - pair orientation
%   I - current
%   WL - width and length
%   T - thickness

Itot=I*T/3.6; %Current multiplied by # of windings
c=1; %proper units for cgs?

X=x-x0;
Y=y-y0;
Z=z-z0;

if n0(1)==1
    B1=((X-W)/sqrt((Z-L)^2+Y^2+(X-W)^2)-(X+W)/sqrt((X+W)^2+Y^2+(Z-L)^2))/sqrt((Z-L)^2+Y^2);
    B2=((Z-L)/sqrt((X-W)^2+Y^2+(Z-L)^2)-(Z+L)/sqrt((X-W)^2+Y^2+(Z+L)^2))/sqrt((X-W)^2+Y^2);
    B3=-((X-W)/sqrt((Z+L)^2+Y^2+(X-W)^2)-(X+W)/sqrt((X+W)^2+Y^2+(Z+L)^2))/sqrt((Z+L)^2+Y^2);
    B4=-((Z-L)/sqrt((X+W)^2+Y^2+(Z-L)^2)-(Z+L)/sqrt((X+W)^2+Y^2+(Z+L)^2))/sqrt((X+W)^2+Y^2);
    
    Bout=-Itot/c*[(B1+B3)*cos(atan(Y/Z))-(B2+B4)*sin(atan(Y/X)),-(B1+B3)*sin(atan(Y/Z))+(B2+B4)*cos(atan(Y/X)),0];
else
    B1=((Y-W)/sqrt((Z-L)^2+X^2+(Y-W)^2)-(Y+W)/sqrt((Y+W)^2+X^2+(Z-L)^2))/sqrt((Z-L)^2+X^2);
    B2=((Z-L)/sqrt((Y-W)^2+X^2+(Z-L)^2)-(Z+L)/sqrt((Y-W)^2+X^2+(Z+L)^2))/sqrt((Y-W)^2+X^2);
    B3=-((Y-W)/sqrt((Z+L)^2+X^2+(Y-W)^2)-(Y+W)/sqrt((Y+W)^2+X^2+(Z+L)^2))/sqrt((Z+L)^2+X^2);
    B4=-((Z-L)/sqrt((Y+W)^2+X^2+(Z-L)^2)-(Z+L)/sqrt((Y+W)^2+X^2+(Z+L)^2))/sqrt((Y+W)^2+X^2);
    
    Bout=-Itot/c*[(B1+B3)*cos(atan(X/Z))-(B2+B4)*sin(atan(X/Y)),-(B1+B3)*sin(atan(X/Z))+(B2+B4)*cos(atan(X/Y)),0];
end


end

