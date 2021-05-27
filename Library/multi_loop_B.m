function Bout = multi_loop_B(x,y,z,x0,y0,z0,n,I,Rin,Rout,H)

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

rp=InvTrans*([dx dy dz]');

%Integral!
C=quadv(@(R)quadv(@(D)B_field_loop( rp(1),rp(2),rp(3),D,R ),-H/2,H/2),Rin,Rout);

%Integral!(which is very slow)
%C=integral(@(R)integral(@(D)B_field_loop( rp(1),rp(2),rp(3),D,R ),-H/2,H/2,'ArrayValued',true),Rin,Rout,'ArrayValued',true);
C=C*I/((Rout-Rin)*H);

Bout=Trans*C';

end