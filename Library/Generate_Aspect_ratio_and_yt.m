function [yt, Aspect]=Generate_Aspect_ratio_and_yt(x10,y10,z10,n10,R1in,R1out,H1,D1,x20,y20,z20,n20,R2in,R2out,H2,D2)


[ytx, Aspectx]=fminbnd(@(y) -Aspectratio_vs_y(y,x10,y10,z10,n10,R1in,R1out,H1,D1,x20,y20,z20,n20,R2in,R2out,H2,D2),y10,y20); 
Aspect=-Aspectx; 
yt=ytx; 