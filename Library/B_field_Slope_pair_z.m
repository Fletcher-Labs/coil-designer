function [ Bout ] = B_field_Slope_pair_z(  x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
delta=1e-4;
dz=delta;

B_pdz=B_field_pair(x,y,z+dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dz=B_field_pair(x,y,z+2*dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdz=B_field_pair(x,y,z-2*dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dz=B_field_pair(x,y,z-2*dz,x0,y0,z0,n0,I,Rin,Rout,H,D);

Bout=(B_m2dz-8*B_mdz+8*B_pdz-B_p2dz)/12/delta;

end

