clc
clear all

syms S11 S22 S33 S23 S13 S12  E11 E22 E33 E12 E13 E23 ...
    F11 F12 F13 F21 F22 F23 F31 F32 F33

ndim = 3; 

S_vect = [S11;S22; S33; S23; S13; S12] ; 
E_vect = [E11;E22; E33; 2*E23; 2*E13; 2*E12] ;
F_vect = [F11;F22; F33; F23; F13; F12 ; F32 ; F31; F21] ;


S_mat = [S11 S12 S13
         S12 S22 S23
         S13 S23 S33] ; 
E_mat = [E11 E12 E13
         E12 E22 E23
         E13 E23 E33] ; 
     
     F_mat = [F11 F12 F13
         F21 F22 F23
         F31 F32 F33] ; 
          
     

GammaS = vectTOmatSYMstress(ndim);


S = sym(zeros(3,3)) ; 

for idim = 1:length(S_vect)
    S = S + squeeze(GammaS(:,:,idim))*S_vect(idim) ;
end

disp('Checking Gamma (correct if all zero)')
S-S_mat

%%%%%%%%%%%%%%%%

Sv = sym(zeros(6,1)) ; 

GammaSinv = matTOvectSYMstress(ndim) ;

for idim = 1:ndim
    for jdim = 1:ndim 
    Sv = Sv + squeeze(GammaSinv(:,idim,jdim))*S_mat(idim,jdim) ;
    end
end

disp('Checking GammaINV (correct if all zero)')
Sv-S_vect


%%%%%%%%%%%%%%%% STRAIN

GammaE = vectTOmatSYMstrain(ndim);


E = sym(zeros(3,3)) ; 

for idim = 1:length(E_vect)
    E = E + squeeze(GammaE(:,:,idim))*E_vect(idim) ;
end

disp('Checking GammaE (correct if all zero)')
E-E_mat

%%%%%%%%%%%%%%%%

Ev = sym(zeros(6,1)) ; 

GammaEinv = matTOvectSYMstrain(ndim) ;

for idim = 1:ndim
    for jdim = 1:ndim 
    Ev = Ev + squeeze(GammaEinv(:,idim,jdim))*E_mat(idim,jdim) ;
    end
end

disp('Checking GammaEINV (correct if all zero)')
Ev-E_vect


%%%%%%%%%%%%%%%% F

Lambda = vectTOmat(ndim);


F = sym(zeros(3,3)) ; 

for idim = 1:length(F_vect)
    F = F + squeeze(Lambda(:,:,idim))*F_vect(idim) ;
end

disp('Checking Lambda (correct if all zero)')
F-F_mat

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

Fv = sym(zeros(9,1)) ; 

LambdaINV = matTOvect(ndim) ;

for idim = 1:ndim
    for jdim = 1:ndim 
    Fv = Fv + squeeze(LambdaINV(:,idim,jdim))*F_mat(idim,jdim) ;
    end
end

disp('Checking GammaEINV (correct if all zero)')
Fv-F_vect
