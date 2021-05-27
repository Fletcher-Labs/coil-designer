function Gout = B_field_Slope_antiHelmholtz_x(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%%%%  Given a pair of anti helmholtz coils at (x0,y0,z0) direct to n0
%%%%  (upper one), with current I, inner radius Rin, outer radius Rout, and
%%%%   thickness H and distance 2D, calculate dBx/dx, dBy/dy, dBz/dz at x,y,z.


%set dx
delta=1e-4;
dx=delta;

%d(Bx)/dx
B_pdx=multi_antiHelmholtz_Bx(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dx=multi_antiHelmholtz_Bx(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdx=multi_antiHelmholtz_Bx(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dx=multi_antiHelmholtz_Bx(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);

Gout=(-B_p2dx+8*B_pdx-8*B_mdx+B_m2dx)/12/delta;

end


