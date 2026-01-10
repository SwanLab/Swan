clc
clear all
%
nstrain =6;
COMPUTE_again = 1 ; 
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
Ivol = sym(zeros(nstrain))  ;
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


NameDiary = 'Ivol_diary.txt' ;
StrInputVariable_symb = 'Cb' ;
StrInputVariable_vect = 'Cb' ;
NameRowsINDEXES = 'ROWS' ;
InputVariableSymbolic = Ivol ;
OutPutVariableVector = 'Ivol'  ;



WriteSymbolicExpressionsMatrix(NameDiary,nstrain,StrInputVariable_symb,StrInputVariable_vect,...
    NameRowsINDEXES,InputVariableSymbolic,OutPutVariableVector)  ;





