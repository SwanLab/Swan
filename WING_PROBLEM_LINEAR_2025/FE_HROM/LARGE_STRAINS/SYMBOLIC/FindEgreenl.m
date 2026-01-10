clc
clear all

syms     F1 F2 F3 F4 F5 F6 F7 F8 F9

 
F_mat = [F1  F6 F5
         F9  F2 F4
        F8   F7 F3 ] ;
ident = [1 0 0; 
         0 1 0 
         0  0  1] ; 

E = 0.5*(F_mat.'*F_mat-ident)

E1 = (E(1,1)) 
E2 = E(2,2)
E3 = E(3,3)
E4 = 2*E(2,3)
E5 = 2*E(1,3)
E6 = 2*E(1,2)