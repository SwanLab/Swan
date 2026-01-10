clc
clear all
%
nstrain =6;
COMPUTE_again = 1;

addpath('../SYMBOLIC/')

if nstrain == 3
    ndim = 2;
    syms  Cb1 Cb2 Cb3
    C_mat = [Cb1  Cb3
        Cb3  Cb2] ;
    
else
    ndim = 3;
    syms  Cb1 Cb2 Cb3 Cb4 Cb5 Cb6
    C_mat = [Cb1 Cb6 Cb5; Cb6 Cb2 Cb4 ; Cb5 Cb4 Cb3] ;
    
end

GammaINV = matTOvectSYMstress(ndim) ;
GammaB = vectTOmatSYMstrain(ndim);
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
    save('Isym_ws.mat','Isym')
else
    load('Isym_ws.mat','Isym')
end



NameDiary = 'Isym_diary.txt' ;
StrInputVariable_symb = 'Cb' ;
StrInputVariable_vect = 'Cb' ;
NameRowsINDEXES = 'ROWS' ;
InputVariableSymbolic = Isym ;
OutPutVariableVector = 'Isym'  ;



WriteSymbolicExpressionsMatrix(NameDiary,nstrain,StrInputVariable_symb,StrInputVariable_vect,...
    NameRowsINDEXES,InputVariableSymbolic,OutPutVariableVector)  ;





