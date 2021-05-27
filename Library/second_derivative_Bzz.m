function ddb = second_derivative_Bzz ( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% set dx
delta=1e-4;
dz=delta;

B_pdz=B_field_pair(x,y,z+dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dz=B_field_pair(x,y,z+2*dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_0=B_field_pair(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdz=B_field_pair(x,y,z-dz,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dz=B_field_pair(x,y,z-2*dz,x0,y0,z0,n0,I,Rin,Rout,H,D);

ddb=(-B_m2dz+16*B_mdz-30*B_0+16*B_pdz-B_p2dz)/12/delta^2;

end

