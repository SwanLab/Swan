clc
clear all
% See 

syms  Finv11 Finv22 Finv33 Finv12 Finv21 Finv13 Finv31 Finv23 Finv32
ndim = 3;
GENERATE= 1;
if ndim == 3   
    Finv= [Finv11 Finv12 Finv13 ;     Finv21 Finv22 Finv23 ;     Finv31 Finv32 Finv33] ;   nstrain  = 6;
    
    
    
    
else
    Finv = [Finv11 Finv12 ;          Finv21 Finv22] ;       nstrain = 3;
   
end

% \begin{equation}
% \label{eq:82eeeeed}
%   \voigt{\e}_c  =   \overbrace{\Par{\matTOvectE{c}{a}{b} [\F^{-1}]_{Aa}  [\F^{-1}]_{Bb}}   \vectTOmatE{A}{B}{C}  }^{\ddSS \R^{eE}} \voigt{\E}_C
% \end{equation}


if GENERATE ==1
    
  %  \matTOvectE{c}{a}{b} [\F^{-1}]_{Aa}  [\F^{-1}]_{Bb}}   \vectTOmatE{A}{B}{C}  

    
    GammaE = vectTOmatSYMstrain(ndim) ;    R = sym(zeros(nstrain,nstrain)) ;
    GammaEinv = matTOvectSYMstrain(ndim) ; 
    for c = 1:nstrain
        for  a= 1:ndim
            for b = 1:ndim 
                for A = 1:ndim 
                    for B = 1:ndim 
                        for C = 1:nstrain
                            R(c,C) = R(c,C) + GammaEinv(c,a,b)*Finv(A,a)*Finv(B,b)*GammaE(A,B,C) ; 
                        end
                    end
                end
            end
        end
    end
    
    save('Rmap.mat','R') ;
    
   
      
else
    load('Rmap.mat')
end




%%%%%

% 
% % Equivalence
% nF = 9 ;
% Fequiv = cell(nF,2) ;
% for i=1:nF
%     Fequiv{i,1} = ['F',num2str(i)];
%     Fequiv{i,2} = ['FgradST(FROWS{',num2str(i),'})'];
% end
% nstrain = 6 ;
% Cequiv = cell(nstrain,nstrain,2) ;
% for istrain =1:nstrain
%     for jstrain = 1:nstrain
%         Cequiv{istrain,jstrain,1} = ['C',num2str(istrain),num2str(jstrain)] ;
%         Cequiv{istrain,jstrain,2} = ['celastST(CROWS{',num2str(istrain),'},',num2str(jstrain),')'] ;
%     end
% end
% 
% diary('TMP_Cmat.txt')
% 
% for i=1:size(CmatFIN,1)
%     for j=i:size(CmatFIN,2)
%         label_i = ['FROWS{',num2str(i),'}'] ;
%         label_j = num2str(j) ;
%         C_loc= ['celasLARGEmat(',label_i,',',label_j,') = '] ;
%         Cstr = sym2str(CmatFIN(i,j)) ;
%         
%         %    if  j >=i
%         for ireplace = 1:size(Fequiv,1)
%             Cstr = strrep(Cstr,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
%         end
%         for ireplace = 1:size(Cequiv,1)
%             for jreplace= 1:size(Cequiv,2)
%                 Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
%             end
%         end
%         
%         %   else
%         %                 Cstr= ['celasLARGEmat(',label_j,',',label_i,') = '] ;
%         
%         %  end
%         disp([C_loc,Cstr,';']) ;
%         
%     end
% end
% 
% for i=1:size(CmatFIN,1)
%     for j=1:size(CmatFIN,2)
%         label_i = ['FROWS{',num2str(i),'}'] ;
%         label_j = num2str(j) ;
%         C_loc= ['celasLARGEmat(',label_i,',',label_j,') = '] ;
%         Cstr = sym2str(CmatFIN(i,j)) ;
%         
%         if  j >=i
%             %         for ireplace = 1:size(Fequiv,1)
%             %            Cstr = strrep(Cstr,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
%             %         end
%             %         for ireplace = 1:size(Cequiv,1)
%             %             for jreplace= 1:size(Cequiv,2)
%             %                 Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
%             %             end
%             %         end
%             
%         else
%             label_i = num2str(i) ; ;
%             label_j =['FROWS{',num2str(j),'}'] ;
%             Cstr= ['celasLARGEmat(',label_j,',',label_i,') '] ;
%             disp([C_loc,Cstr,';']) ;
%         end
%         
%         
%     end
% end
% 
% 
% diary off
% 
% 
