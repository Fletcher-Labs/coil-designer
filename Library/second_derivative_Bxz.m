function [ Bout ] = second_derivative_Bxz( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
delta=1e-4;
dx=delta;

B_pdx=B_field_Slope_pair_z(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dx=B_field_Slope_pair_z(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdx=B_field_Slope_pair_z(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dx=B_field_Slope_pair_z(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);

Bout=(B_m2dx-8*B_mdx+8*B_pdx-B_p2dx)/12/delta;

end

