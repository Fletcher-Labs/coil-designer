function Bout = B_field_loop( xp,yp,zp,D,R )
%%%%  This function should be used to calculate an ideal current loop with identity current which is located in Z=D and direct to [0,0,1] %%%%
%%%% The output is the B-field vector at xp,yp,zp %%%%
%%%% units in cm and gauss and amps  %%%%

Mu=4*pi()/10;

zp=zp-D;
%Convert to cylindrical coordinate system
rho=sqrt(xp^2+yp^2);    
fai=atan2(yp,xp);         

%计算磁场的几何因子
k=((4.*R.*rho)./((R+rho).^2+zp.^2)).^(1/2);%椭圆积分函数的模数
[K,E]=ellipke(k.^2);            %计算第一类，第二类完全椭圆积分K,E

if rho==0
    Crho=0.*R;   %when rho=0, B_rho is 0
else
    Crho=(zp./(rho.*((R+rho).^2+zp.^2).^(1/2))).*(-K+E.*(R.^2+rho.^2+zp.^2)./((R-rho).^2+zp.^2)); 
end;

Bx=Crho.*cos(fai).*(Mu/(2*pi())); 
By=Crho.*sin(fai).*(Mu/(2*pi())); 
Bz=(1./sqrt((R+rho).^2+zp.^2)).*(K+E.*(R.^2-rho.^2-zp.^2)./((R-rho).^2+zp.^2)).*(Mu/(2*pi())); 

Bout=[Bx By Bz];
    
end

