function [ Bout ] = Rectangular_B_field_pair( X,Y,Z,x0,y0,z0,n0,I,W,L,T,S )
% Calculates the magnetic field due to a rectangular pair of coils
%   XYZ - point of calculation
%   x0y0z0 - coil pair center
%   n0 - pair orientation
%   I - current
%   WL - width and length (half)
%   TS - thickness (full) and separation (half)

Bu=Rectangular_B_field(X,Y,Z,x0+S*n0(1),y0+S*n0(2),z0+S*n0(3),n0,I,W,L,T); %field of upper coil
Bd=Rectangular_B_field(X,Y,Z,x0-S*n0(1),y0-S*n0(2),z0-S*n0(3),n0,I,W,L,T); %field of lower coil
Bout=Bu+Bd; 

end

