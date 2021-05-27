function [yt, AR]=Transitionpoint_Aspect_ratio_and_yt(x10,y10,z10,n10,R1in,R1out,H1,D1,x20,y20,z20,n20,R2in,R2out,H2,D2)

%%%%  calculate the transition point determined by the dimension and
%%%%  geometry, and give the aspect ration of the transition point

%%%% Pair1: centered at (X10,Y10,Z10), direct to n10 with R1in, R1out, H1
%%%% and D1
%%%% Pair2: centered at (X20,Y20,Z20), direct to n20 with R2in, R2out, H2
%%%% and D2

[ym, A]=fminbnd(@(y) -Aspectratio_vs_y(y,x10,y10,z10,n10,R1in,R1out,H1,D1,x20,y20,z20,n20,R2in,R2out,H2,D2),y10,y20); %使用matlab中的fminbnd函数求极小值
AR = -A;
yt=ym; 

