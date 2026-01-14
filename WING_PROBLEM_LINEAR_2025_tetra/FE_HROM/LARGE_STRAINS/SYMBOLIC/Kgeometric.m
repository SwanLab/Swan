clc
clear all

syms  S1 S2 S3 S4 S5 S6
ndim = 3;  
S_mat = [S1 S6 S5
         S6 S2 S4
         S5 S4 S3] ; 
     
     if ndim == 2
S_mat = [S1 S3 
         S3 S2] ;      
     end
    
 GammaS = vectTOmatSYMstress(ndim) ;
 Lambda = vectTOmat(ndim);
 LambdaINV = matTOvect(ndim) ;
 
 celasGEO = sym(zeros(ndim^2))  ; 
 nmatnon = ndim^2 ;
 nstrain  = 6  ; 
 if ndim ==2
     nstrain =3 ; 
 end
 
 for   a = 1:nmatnon
     for  e = 1:nmatnon
          disp(['C(',num2str(a),',',num2str(e),')'])
         for c= 1:ndim
                 for C = 1:ndim
                     for B= 1:ndim
                       celasGEO(a,e) =  celasGEO(a,e)  +    LambdaINV(a,c,B)*Lambda(c,C,e)*S_mat(C,B) ;
                     end
                 end             
         end
     end
 end
 
 save('celasGEO.mat','celasGEO')

     
% 
% 
% S = sym(zeros(3,3)) ; 
% 
% for idim = 1:length(S_vect)
%     S = S + squeeze(GammaS(:,:,idim))*S_vect(idim) ;
% end
% 
% disp('Checking Gamma (correct if all zero)')
% S-S_mat
% 
% %%%%%%%%%%%%%%%%
% 
% Sv = sym(zeros(6,1)) ; 
% 
% GammaSinv = matTOvectSYMstress(ndim) ;
% 
% for idim = 1:ndim
%     for jdim = 1:ndim 
%     Sv = Sv + squeeze(GammaSinv(:,idim,jdim))*S_mat(idim,jdim) ;
%     end
% end
% 
% disp('Checking GammaINV (correct if all zero)')
% Sv-S_vect
% 
% 
% %%%%%%%%%%%%%%%% STRAIN
% 
% GammaE = vectTOmatSYMstrain(ndim);
% 
% 
% E = sym(zeros(3,3)) ; 
% 
% for idim = 1:length(E_vect)
%     E = E + squeeze(GammaE(:,:,idim))*E_vect(idim) ;
% end
% 
% disp('Checking GammaE (correct if all zero)')
% E-E_mat
% 
% %%%%%%%%%%%%%%%%
% 
% Ev = sym(zeros(6,1)) ; 
% 
% GammaEinv = matTOvectSYMstrain(ndim) ;
% 
% for idim = 1:ndim
%     for jdim = 1:ndim 
%     Ev = Ev + squeeze(GammaEinv(:,idim,jdim))*E_mat(idim,jdim) ;
%     end
% end
% 
% disp('Checking GammaEINV (correct if all zero)')
% Ev-E_vect
% 
% 
% %%%%%%%%%%%%%%%% F
% 
% Lambda = vectTOmat(ndim);
% 
% 
% F = sym(zeros(3,3)) ; 
% 
% for idim = 1:length(F_vect)
%     F = F + squeeze(Lambda(:,:,idim))*F_vect(idim) ;
% end
% 
% disp('Checking Lambda (correct if all zero)')
% F-F_mat
% 
% %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%
% 
% Fv = sym(zeros(9,1)) ; 
% 
% LambdaINV = matTOvect(ndim) ;
% 
% for idim = 1:ndim
%     for jdim = 1:ndim 
%     Fv = Fv + squeeze(LambdaINV(:,idim,jdim))*F_mat(idim,jdim) ;
%     end
% end
% 
% disp('Checking GammaEINV (correct if all zero)')
% Fv-F_vect
