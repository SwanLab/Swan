clc
clear all

%syms  C1 C2 C3 C4 C5 C6
COMPUTE_again = 1; 
ndim = 3;
C_mat = eye(3) ;

if ndim == 2
    C_mat = [C1 C3
        C3 C2] ;
end

GammaINV = matTOvectSYMstress(ndim) ;
Gamma = vectTOmatSYMstress(ndim) ;

GammaB = vectTOmatSYMstrain(ndim)
 
nstrain  = 6  ;
if ndim ==2
    nstrain =3 ;
end
Ivol = sym(zeros(nstrain))  ;

% G = sym(zeros(nstrain))  ; 
% 
% for   c = 1:nstrain
%     for  b = 1:nstrain
%         disp(['G(',num2str(c),',',num2str(b),')'])
%         for C= 1:ndim
%             for D =  1:ndim                
%                 
%                     G(c,b) =  G(c,b)  +  Gamma(C,D,c)*GammaB(C,D,b) ;   
%                      
%             end
%         end
%     end
% end

if COMPUTE_again == 1 
for   a = 1:nstrain
    for  b = 1:nstrain
        disp(['Ivol(',num2str(a),',',num2str(b),')'])
        for A= 1:ndim
            for B = 1:ndim                
                for C= 1:ndim
                    for D = 1:ndim 
                    Ivol(a,b) =  Ivol(a,b)  +  GammaINV(a,A,B)*C_mat(A,B)*C_mat(C,D)*GammaB(C,D,b) ;   
                    end
                end
            end
        end
    end
end

save('Ivol_ws.mat','Ivol')

else
   load('Ivol_ws.mat','Ivol') 
end

