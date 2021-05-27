function Gout = B_field_Slope_loop_z( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H )
%%%%  Given a coil at (x0,y0,z0) direct to n0, with current I, inner radius
%%%%  Rin, outer radius Rout, and thickness H, calculate dBx/dx, dBy/dy,
%%%%  dBz/dz at x,y,z.

%set dx
delta=1e-4;
dx=delta;

%d(Bz)/dz
B_pdz=multi_loop_B_z(x,y,z+dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dz=multi_loop_B_z(x,y,z+2*dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdz=multi_loop_B_z(x,y,z-dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dz=multi_loop_B_z(x,y,z-2*dx,x0,y0,z0,n0,I,Rin,Rout,H);

Gout=(-B_p2dz+8*B_pdz-8*B_mdz+B_m2dz)/12/delta;

end

