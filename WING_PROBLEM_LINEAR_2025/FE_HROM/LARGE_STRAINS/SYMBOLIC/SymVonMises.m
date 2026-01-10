clc
clear all

 
ndim = 3;
if ndim == 3
    
    
    
    
    
else
    nF = 4; nstrain =3 ; 
    F11 = F1;  F22 = F2 ;  F12 = F3 ; F21 = F4 ;
    F = [F11 F12 ;          F21 F22] ;
    
    %
    P11 = P1;  P22 = P2 ;  P12 = P3 ; P21 = P4 ;
    P = [P11 P12 ;          P21 P22] ;
    
    
    sigmaS = F.'*P/detF ;
    
    sigma{1} = sigmaS(1,1) ;
    sigma{2}= sigmaS(2,2) ;
    sigma{3} = sigmaS(1,2) ;
end





%%%%%


% Equivalence

Fequiv = cell(nF,2) ;
for i=1:nF
    Fequiv{i,1} = ['F',num2str(i)];
    Fequiv{i,2} = ['FgradST(FROWS{',num2str(i),'})'];
end

Pequiv = cell(nF,2) ;
for i=1:nF
    Pequiv{i,1} = ['P',num2str(i)];
    Pequiv{i,2} = ['PoneST(FROWS{',num2str(i),'})'];
end



% CauchyEquiv = cell(nstrain,2) ;
% for istrain =1:nstrain
%     CauchyEquiv{istrain,1} = ['sigma',num2str(istrain)] ;
%     CauchyEquiv{istrain,2} = ['CauchyStress(CROWS{',num2str(istrain),'})'] ;
% end

diary('CuachyStress.txt')


for istrain = 1:nstrain
    
    
    label_i = ['CROWS{',num2str(istrain),'}'] ;
    C_loc= ['CauchyStress(',label_i,') = '] ;
    Cstr = sym2str(sigma{istrain}) ;
    
    for ireplace = 1:size(Fequiv,1)
        Cstr = strrep(Cstr,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
    end
    
    for ireplace = 1:size(Fequiv,1)
        Cstr = strrep(Cstr,Pequiv{ireplace,1},Pequiv{ireplace,2}) ;
    end
    
    
    disp([C_loc,Cstr,';']) ;
    
end







%
% detF = sym2str(detF) ;
% for ireplace = 1:size(Fequiv,1)
%     detF = strrep(detF,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
% end
% %         for ireplace = 1:size(Cequiv,1)
% %             for jreplace= 1:size(Cequiv,2)
% %                 Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
% %             end
% %         end
%
% %   else
% %                 Cstr= ['celasLARGEmat(',label_j,',',label_i,') = '] ;
%
% %  end
% disp(['detF =',detF,';']) ;
%
% %     end
% % end

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


diary off


