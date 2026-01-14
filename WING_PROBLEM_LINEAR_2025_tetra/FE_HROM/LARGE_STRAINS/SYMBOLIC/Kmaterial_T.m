clc
clear all

syms  F1 F2 F3 F4 F5 F6 F7 F8 F9
ndim = 3;
GENERATE= 1;
if ndim == 3
    F11 = F1; F22 = F2 ; F33 = F3;  F23 = F4; F13 = F5; F12 = F6 ; F32 = F7; F31 = F8 ; F21 =F9 ;
    F = [F11 F12 F13 ;     F21 F22 F23 ;     F31 F32 F33] ;   nstrain  = 6;
    
    syms C11 C12 C13 C14 C15 C16
    syms C21 C22 C23 C24 C25 C26
    syms C31 C32 C33 C34 C35 C36
    syms C41 C42 C43 C44 C45 C46
    syms C51 C52 C53 C54 C55 C56
    syms C61 C62 C63 C64 C65 C66
    
    Celas = [C11 C12 C13 C14 C15 C16
        C21 C22 C23 C24 C25 C26
        C31 C32 C33 C34 C35 C36
        C41 C42 C43 C44 C45 C46
        C51 C52 C53 C54 C55 C56
        C61 C62 C63 C64 C65 C66] ;
    
    
    
else
    F11 = F1;  F22 = F2 ;  F12 = F3 ; F21 = F4 ;
    F = [F11 F12 ;          F21 F22] ;       nstrain = 3;
    syms C11 C12 C13
    syms C21 C22 C23
    syms C31 C32 C33
    Celas=[C11,C12,C13
        C21 C22 C23
        C31 C32 C33] ;
end

if GENERATE ==1
    GammaE = vectTOmatSYMstrain(ndim) ;  GammaS = vectTOmatSYMstress(ndim) ;
    Lambda = vectTOmat(ndim); LambdaINV = matTOvect(ndim) ; T = sym(zeros(ndim^2,nstrain)) ;
    for a = 1:ndim^2
        for  h= 1:nstrain
            disp(['T(',num2str(a),',',num2str(h),')'])
            for b = 1:ndim
                for B = 1:ndim
                    for D = 1:ndim
                        T(a,h) = T(a,h) + LambdaINV(a,b,B)*F(b,D)*GammaS(D,B,h) ;
                    end
                end
            end
        end
    end
    
    save('TransF.mat','T') ;
    
    GammaEinv = matTOvectSYMstrain(ndim) ;
    Tbar = sym(zeros(nstrain,ndim^2)) ;
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
    
    T-Tbar.'
    
    
    
    
    CmatFIN = T*Celas*T.'
    save('TMP_Cmat','CmatFIN')
    
else
    load('TMP_Cmat')
end




%%%%%


% Equivalence
nF = 9 ;
Fequiv = cell(nF,2) ;
for i=1:nF
    Fequiv{i,1} = ['F',num2str(i)];
    Fequiv{i,2} = ['FgradST(FROWS{',num2str(i),'})'];
end
nstrain = 6 ;
Cequiv = cell(nstrain,nstrain,2) ;
for istrain =1:nstrain
    for jstrain = 1:nstrain
        Cequiv{istrain,jstrain,1} = ['C',num2str(istrain),num2str(jstrain)] ;
        Cequiv{istrain,jstrain,2} = ['celastST(CROWS{',num2str(istrain),'},',num2str(jstrain),')'] ;
    end
end

diary('TMP_Cmat.txt')

for i=1:size(CmatFIN,1)
    for j=i:size(CmatFIN,2)
        label_i = ['FROWS{',num2str(i),'}'] ;
        label_j = num2str(j) ;
        C_loc= ['celasLARGEmat(',label_i,',',label_j,') = '] ;
        Cstr = sym2str(CmatFIN(i,j)) ;
        
        %    if  j >=i
        for ireplace = 1:size(Fequiv,1)
            Cstr = strrep(Cstr,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
        end
        for ireplace = 1:size(Cequiv,1)
            for jreplace= 1:size(Cequiv,2)
                Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
            end
        end
        
        %   else
        %                 Cstr= ['celasLARGEmat(',label_j,',',label_i,') = '] ;
        
        %  end
        disp([C_loc,Cstr,';']) ;
        
    end
end

for i=1:size(CmatFIN,1)
    for j=1:size(CmatFIN,2)
        label_i = ['FROWS{',num2str(i),'}'] ;
        label_j = num2str(j) ;
        C_loc= ['celasLARGEmat(',label_i,',',label_j,') = '] ;
        Cstr = sym2str(CmatFIN(i,j)) ;
        
        if  j >=i
            %         for ireplace = 1:size(Fequiv,1)
            %            Cstr = strrep(Cstr,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
            %         end
            %         for ireplace = 1:size(Cequiv,1)
            %             for jreplace= 1:size(Cequiv,2)
            %                 Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
            %             end
            %         end
            
        else
            label_i = num2str(i) ; ;
            label_j =['FROWS{',num2str(j),'}'] ;
            Cstr= ['celasLARGEmat(',label_j,',',label_i,') '] ;
            disp([C_loc,Cstr,';']) ;
        end
        
        
    end
end


diary off


