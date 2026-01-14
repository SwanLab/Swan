clc
clear all

nstrain = 6

if nstrain ==3
    syms E1 E2 E3
    ndim = 2;
    
    E =  sym(zeros(ndim,ndim)) ;
    
    E(1,1) = E1;
    E(1,2) = 0.5*E3;
    E(2,1) = E(1,2) ;
    E(2,2)  = E2 ;
    
    C = 2*E + sym(eye(ndim)) ;
    invC = inv(C) ;
    
    detC =det(C) ;
    
    invCdet = sym(zeros(nstrain,1)) ;
    
    invCdet(1) = invC(1,1)*detC
    invCdet(2) = invC(2,2)*detC
    invCdet(3) = invC(1,2)*detC
    
    
    %      E_equiv = cell(nstrain,2) ;
    %     for i=1:nstrain
    %         E_equiv{i,1} = ['E',num2str(i)];
    %         E_equiv{i,2} = ['E(FROWS{',num2str(i),'})'];
    %     end
    
    
else
    
    
    syms E1 E2 E3 E4 E5 E6 
    ndim = 3;
    
    E =  sym(zeros(ndim,ndim)) ;
    
    E(1,1) = E1;
    E(2,2)  = E2 ;
    E(3,3)  = E3 ;         
    E(2,3) = 0.5*E4;
    E(1,3) = 0.5*E5;
    E(1,2) = 0.5*E6;
    
    E(2,1) = E(1,2) ;
    E(3,1) = E(1,3) ;
    E(3,2) = E(2,3) ;
    
    C = 2*E + sym(eye(ndim)) ;
    invC = inv(C) ;
    
    detC =det(C) ;
    
    invCdet = sym(zeros(nstrain,1)) ;
    
    invCdet(1) = invC(1,1)*detC; 
    invCdet(2) = invC(2,2)*detC;
    invCdet(3) = invC(3,3)*detC;
    invCdet(4) = invC(2,3)*detC;
    invCdet(5) = invC(1,3)*detC;
    invCdet(6) = invC(1,2)*detC;
    
    
    %% 
     
    Eequiv = cell(nstrain,2) ;
    for i=1:nstrain
        Eequiv{i,1} = ['E',num2str(i)];
        Eequiv{i,2} = ['E(SROWS{',num2str(i),'})'];
    end
    
    diary('CrightINVtext.txt')
    
    for istrain = 1:nstrain 
        
        STRLOC = sym2str(invCdet(istrain)) ;
        
        for jstrain = 1:nstrain 
         STRLOC = strrep(STRLOC,Eequiv{jstrain,1},Eequiv{jstrain,2}) ; 
        end
         
        Cstr = ['Cinv(','SROWS{',num2str(istrain),'}',') = (',STRLOC,')./detC; '] ; 
        disp(Cstr)
        
    end
    
    diary off
     
    
end
