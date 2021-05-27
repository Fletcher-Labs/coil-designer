function ddb = second_derivative_norm_Bxx ( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% set dx
delta=1e-4;
dx=delta;

B_pdx=norm(B_field_pair(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D));
B_p2dx=norm(B_field_pair(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D));
B_0=norm(B_field_pair(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D));
B_mdx=norm(B_field_pair(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D));
B_m2dx=norm(B_field_pair(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D));

ddb=(-B_m2dx+16*B_mdx-30*B_0+16*B_pdx-B_p2dx)/12/delta^2;

end
