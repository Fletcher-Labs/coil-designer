function Gout = B_field_Slope_loop_y( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H )
%%%%  Given a coil at (x0,y0,z0) direct to n0, with current I, inner radius
%%%%  Rin, outer radius Rout, and thickness H, calculate dBx/dx, dBy/dy,
%%%%  dBz/dz at x,y,z.

%set dx
delta=1e-4;
dx=delta;

%d(By)/dy
B_pdy=multi_loop_B_y(x,y+dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dy=multi_loop_B_y(x,y+2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdy=multi_loop_B_y(x,y-dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dy=multi_loop_B_y(x,y-2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H);

Gout=(-B_p2dy+8*B_pdy-8*B_mdy+B_m2dy)/12/delta;



end

