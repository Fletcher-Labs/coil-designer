function Bout = multi_loop_B_x(x,y,z,x0,y0,z0,n,I,Rin,Rout,H)

%%%%  Calculate the magnetic field of a coil whose center is located in (x0,y0,z0) and direct to n, and have a current of I %%%%

n=n./norm(n);

if (abs(n(1))==1)   %find a new vector which is vertical to n
    l=[n(3),0,-n(1)];
else
    l=[0,n(3),-n(2)];
end
l=l/norm(l);
m=cross(n,l);            %m is the cross product of n and l
Trans=[l',m',n'];           %Transformation matrix from coil coordinator to lab coordinator
InvTrans=inv(Trans);       

dx=x-x0;
dy=y-y0;
dz=z-z0;

xp=dx*InvTrans(1,1)+dy*InvTrans(1,2)+dz*InvTrans(1,3); 
yp=dx*InvTrans(2,1)+dy*InvTrans(2,2)+dz*InvTrans(2,3);
zp=dx*InvTrans(3,1)+dy*InvTrans(3,2)+dz*InvTrans(3,3);

%Integral!
Cx=quad2d(@(D,R) B_field_loop_x(xp,yp,zp,D,R),-H/2,H/2,Rin,Rout);
Cy=quad2d(@(D,R) B_field_loop_y(xp,yp,zp,D,R),-H/2,H/2,Rin,Rout);
Cz=quad2d(@(D,R) B_field_loop_z(xp,yp,zp,D,R),-H/2,H/2,Rin,Rout);

Bx=(Cx*Trans(1,1)+Cy*Trans(1,2)+Cz*Trans(1,3));

Bout=Bx*I/((Rout-Rin)*H);

end