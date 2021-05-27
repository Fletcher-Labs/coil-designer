function Gout = B_field_Slope_antiHelmholtz_z(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%%%%  Given a pair of anti helmholtz coils at (x0,y0,z0) direct to n0
%%%%  (upper one), with current I, inner radius Rin, outer radius Rout, and
%%%%   thickness H and distance 2D, calculate dBx/dx, dBy/dy, dBz/dz at x,y,z.


%set dx
delta=1e-4;
dx=delta;

%d(Bz)/dz
B_pdz=multi_antiHelmholtz_Bz(x,y,z+dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dz=multi_antiHelmholtz_Bz(x,y,z+2*dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdz=multi_antiHelmholtz_Bz(x,y,z-dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dz=multi_antiHelmholtz_Bz(x,y,z-2*dx,x0,y0,z0,n0,I,Rin,Rout,H,D);

Gout=(-B_p2dz+8*B_pdz-8*B_mdz+B_m2dz)/12/delta;

end