function Gout = B_field_Slope_antiHelmholtz(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%%%%  Given a pair of anti helmholtz coils at (x0,y0,z0) direct to n0
%%%%  (upper one), with current I, inner radius Rin, outer radius Rout, and
%%%%   thickness H and distance 2D, calculate dBx/dx, dBy/dy, dBz/dz at x,y,z.


%set dx
delta=1e-4;
dx=delta;

%d(Bx)/dx
B_pdx=multi_antiHelmholtz_B(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dx=multi_antiHelmholtz_B(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdx=multi_antiHelmholtz_B(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dx=multi_antiHelmholtz_B(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D);

slopeBx=(-B_p2dx+8*B_pdx-8*B_mdx+B_m2dx)/12/delta;
slopeBx=slopeBx(1);

%d(By)/dy
B_pdy=multi_antiHelmholtz_B(x,y+dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dy=multi_antiHelmholtz_B(x,y+2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdy=multi_antiHelmholtz_B(x,y-dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dy=multi_antiHelmholtz_B(x,y-2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);

slopeBy=(-B_p2dy+8*B_pdy-8*B_mdy+B_m2dy)/12/delta;
slopeBy=slopeBy(2);

%d(Bz)/dz
B_pdz=multi_antiHelmholtz_B(x,y,z+dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dz=multi_antiHelmholtz_B(x,y,z+2*dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdz=multi_antiHelmholtz_B(x,y,z-dx,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dz=multi_antiHelmholtz_B(x,y,z-2*dx,x0,y0,z0,n0,I,Rin,Rout,H,D);

slopeBz=(-B_p2dz+8*B_pdz-8*B_mdz+B_m2dz)/12/delta;
slopeBz=slopeBz(3);

Gout=[slopeBx slopeBy slopeBz];

end


