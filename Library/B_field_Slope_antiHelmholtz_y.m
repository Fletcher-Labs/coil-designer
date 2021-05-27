function Gout = B_field_Slope_antiHelmholtz_y(x,y,z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%%%%  Given a pair of anti helmholtz coils at (x0,y0,z0) direct to n0
%%%%  (upper one), with current I, inner radius Rin, outer radius Rout, and
%%%%   thickness H and distance 2D, calculate dBx/dx, dBy/dy, dBz/dz at x,y,z.


%set dx
delta=1e-4;
dx=delta;
%d(By)/dy
B_pdy=multi_antiHelmholtz_By(x,y+dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_p2dy=multi_antiHelmholtz_By(x,y+2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_mdy=multi_antiHelmholtz_By(x,y-dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);
B_m2dy=multi_antiHelmholtz_By(x,y-2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H,D);

Gout=(-B_p2dy+8*B_pdy-8*B_mdy+B_m2dy)/12/delta;

end


