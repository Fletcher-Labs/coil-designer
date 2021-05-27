function Aspect=Aspectratio_vs_y(y,x10,y10,z10,n10,R1in,R1out,H1,D1,x20,y20,z20,n20,R2in,R2out,H2,D2)

B1y=multi_antiHelmholtz_By(0,y,0,x10,y10,z10,n10,1,R1in,R1out,H1,D1);
B2y=multi_antiHelmholtz_By(0,y,0,x20,y20,z20,n20,1,R2in,R2out,H2,D2);

slope1B=B_field_Slope_antiHelmholtz_z(0,y,0,x10,y10,z10,n10,1,R1in,R1out,H1,D1);
slope2B=B_field_Slope_antiHelmholtz_z(0,y,0,x20,y20,z20,n20,1,R2in,R2out,H2,D2);

%计算x,y方向磁场梯度d(Bx)/dx,d(By)/dy
slope1By=B_field_Slope_antiHelmholtz_y(0,y,0,x10,y10,z10,n10,1,R1in,R1out,H1,D1);
slope2By=B_field_Slope_antiHelmholtz_y(0,y,0,x20,y20,z20,n20,1,R2in,R2out,H2,D2);


%生成系数矩阵
M=[B1y,B2y;slope1B,slope2B];
G=[0,1]';
I=M\G;%计算电流
Aspect=-1/([slope1By,slope2By]*I)-1;%计算展宽比Aspect=Gx/Gy
