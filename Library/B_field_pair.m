function [ Bout ] = B_field_pair( X,Y,Z,x0,y0,z0,n0,I,Rin,Rout,H,D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Bu=multi_loop_B(X,Y,Z,x0+D*n0(1),y0+D*n0(2),z0+D*n0(3),n0,I,Rin,Rout,H); %field of upper coil
Bd=multi_loop_B(X,Y,Z,x0-D*n0(1),y0-D*n0(2),z0-D*n0(3),n0,I,Rin,Rout,H); %field of lower coil
Bout=Bu+Bd; 

end

