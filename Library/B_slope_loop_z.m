function G =B_slope_loop_z( xp,yp,zp,D,R )
%set dx
delta=1e-4;
dx=delta;


B_pdz=B_field_loop_z( xp,yp,zp+dx,D,R );
B_p2dz=B_field_loop_z( xp,yp,zp+2*dx,D,R );
B_mdz=B_field_loop_z( xp,yp,zp-dx,D,R );
B_m2dz=B_field_loop_z( xp,yp,zp-2*dx,D,R );

G=(-B_p2dz+8*B_pdz-8*B_mdz+B_m2dz)/12/delta;
end

