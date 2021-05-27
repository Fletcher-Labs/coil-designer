function Gout = B_field_Slope_loop( x,y,z,x0,y0,z0,n0,I,Rin,Rout,H )
%%%%  Given a coil at (x0,y0,z0) direct to n0, with current I, inner radius
%%%%  Rin, outer radius Rout, and thickness H, calculate dBx/dx, dBy/dy,

%set dx
delta=1e-4;
dx=delta;

%d(Bx)/dx
B_pdx=multi_loop_B(x+dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dx=multi_loop_B(x+2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdx=multi_loop_B(x-dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dx=multi_loop_B(x-2*dx,y,z,x0,y0,z0,n0,I,Rin,Rout,H);

slopeBx=(-B_p2dx(1)+8*B_pdx(1)-8*B_mdx(1)+B_m2dx(1))/12/delta;

%d(By)/dy
B_pdy=multi_loop_B(x,y+dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dy=multi_loop_B(x,y+2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdy=multi_loop_B(x,y-dx,z,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dy=multi_loop_B(x,y-2*dx,z,x0,y0,z0,n0,I,Rin,Rout,H);

slopeBy=(-B_p2dy(2)+8*B_pdy(2)-8*B_mdy(2)+B_m2dy(2))/12/delta;

%d(Bz)/dz
B_pdz=multi_loop_B(x,y,z+dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_p2dz=multi_loop_B(x,y,z+2*dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_mdz=multi_loop_B(x,y,z-dx,x0,y0,z0,n0,I,Rin,Rout,H);
B_m2dz=multi_loop_B(x,y,z-2*dx,x0,y0,z0,n0,I,Rin,Rout,H);

slopeBz=(-B_p2dz(3)+8*B_pdz(3)-8*B_mdz(3)+B_m2dz(3))/12/delta;

Gout=[slopeBx slopeBy slopeBz];

end

