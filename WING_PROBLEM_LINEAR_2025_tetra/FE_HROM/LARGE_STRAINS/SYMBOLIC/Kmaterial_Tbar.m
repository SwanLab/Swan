clc
clear all

syms  F1 F2 F3 F4 F5 F6 F7 F8 F9
ndim = 3;
if ndim == 3
F11 = F1; F22 = F2 ; F33 = F3;  F23 = F4; F13 = F5; F12 = F6 ; F32 = F7; F31 = F8 ; F21 =F9 ;
F = [F11 F12 F13 ;     F21 F22 F23 ;     F31 F32 F33] ;   nstrain  = 6;
else
    F11 = F1;  F22 = F2 ;  F12 = F3 ; F21 = F4 ;  
    F = [F11 F12 ;              F21 F22] ;       nstrain = 3;
end
GammaEinv = matTOvectSYMstrain(ndim) ; 
Lambda = vectTOmat(ndim);   Tbar = sym(zeros(nstrain,ndim^2)) ;
for j = 1:nstrain
    for  e= 1:ndim^2
        disp(['Tbar(',num2str(j),',',num2str(e),')'])
        for C = 1:ndim
            for H = 1:ndim
                for c = 1:ndim        
                    Tbar(j,e) = Tbar(j,e) + GammaEinv(j,C,H)*F(c,H)*Lambda(c,C,e) ;
                end
            end
        end
    end
end

save('TransF_bar.mat','Tbar') ;


 