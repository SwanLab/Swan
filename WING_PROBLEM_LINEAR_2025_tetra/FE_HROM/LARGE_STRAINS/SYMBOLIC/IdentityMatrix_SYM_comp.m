clc
clear all

COMPUTE_again = 1;
ndim = 3;
C_mat =  eye(3);
% if ndim == 2
%     C_mat = [C1 C3
%         C3 C2] ;
% end
 
GammaINV = matTOvectSYMstress(ndim) ;
Gamma = vectTOmatSYMstress(ndim) ;

GammaB = vectTOmatSYMstrain(ndim);

nstrain  = 6  ;
% if ndim ==2
%     nstrain =3 ;
% end
Isym = sym(zeros(nstrain))  ;

 

if COMPUTE_again == 1
    for   a = 1:nstrain
        for  b = 1:nstrain
            disp(['Isym(',num2str(a),',',num2str(b),')'])
            for A= 1:ndim
                for B = 1:ndim
                    for C= 1:ndim
                        for D = 1:ndim
                            Isym(a,b) =  Isym(a,b)  +  (GammaINV(a,A,B)*(C_mat(A,C)*C_mat(B,D)+ C_mat(A,D)*C_mat(B,C))*GammaB(C,D,b))/2 ;
                        end
                    end
                end
            end
        end
    end
    
    save('Isym_ws_comp.mat','Isym')
    
else
    load('Isym_ws_comp.mat','Isym')
end

