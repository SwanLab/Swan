clc
clear all

syms J11 J12 J13 J21 J22 J23 J31 J32 J33

J = [J11 J12 J13; J21 J22 J23; J31 J32 J33]

J = [J11 J12; J21 J22]

invJt = inv(J.')

detJ =det(J) ; 

invJt_detJ = invJt*detJ ; 
simple(invJt_detJ)



