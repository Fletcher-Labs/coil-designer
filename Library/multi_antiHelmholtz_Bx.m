function Bout = multi_antiHelmholtz_Bx( X,Y,Z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%%%%  Calculate the magnetic field in (XYZ) of a pair of anti coils which
%%%%  locate in (x0,y0,z0) with current I, inner radius Rin, Outer radius
%%%%  Rout and thickness H and distance 2D.
%%%%  the upper coil direct to n0.

Bu=multi_loop_B_x(X,Y,Z,x0+D*n0(1),y0+D*n0(2),z0+D*n0(3),n0,I,Rin,Rout,H); %field of upper coil
Bd=multi_loop_B_x(X,Y,Z,x0-D*n0(1),y0-D*n0(2),z0-D*n0(3),-n0,I,Rin,Rout,H); %field of lower coil
Bout=Bu+Bd; 


end

