clc
clear all
% See 

syms  phi
ndim = 2;
nstrain  =3 ; 
Q = [cos(phi) -sin(phi); sin(phi) cos(phi)] ; 
 
     
GammaE = vectTOmatSYMstrain(ndim) ;    R = sym(zeros(nstrain,nstrain)) ;
GammaEinv = matTOvectSYMstrain(ndim) ;
for c = 1:nstrain
    for  a= 1:ndim
        for b = 1:ndim
            for A = 1:ndim
                for B = 1:ndim
                    for C = 1:nstrain
                        R(c,C) = R(c,C) + GammaEinv(c,a,b)*Q(A,a)*Q(B,b)*GammaE(A,B,C) ;
                    end
                end
            end
        end
    end
end