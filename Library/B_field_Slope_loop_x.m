function Gout = B_field_Slope_loop_x( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H )
%%%%  Given a coil at (x0,y0,z0) direct to n0, with current I, inner radius
%%%%  Rin, outer radius Rout, and thickness H, calculate dBx/dx, dBy/dy,
%%%%  dBz/dz at x,y,z.

%set dx
delta=1e-4;
dx=delta;

%d(Bx)/dx
B_pdx=multi_loop_B_x(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dx=multi_loop_B_x(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdx=multi_loop_B_x(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dx=multi_loop_B_x(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);

Gout=(-B_p2dx+8*B_pdx-8*B_mdx+B_m2dx)/12/delta;
end

