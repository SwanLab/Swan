clc
clear all

syms  F1 F2 F3 F4 F5 F6 F7 F8 F9  
%syms  P1 P2 P3 P4 P5 P6 P7 P8 P9 
ndim = 2;
if ndim == 3
    F11 = F1; F22 = F2 ; F33 = F3;  F23 = F4; F13 = F5; F12 = F6 ; F32 = F7; F31 = F8 ; F21 =F9 ;
    F = [F11 F12 F13 ;     F21 F22 F23 ;     F31 F32 F33] ;   
    
    
%     P11 = P1; P22 = P2 ; P33 = P3;  P23 = P4; P13 = P5; P12 = P6 ; P32 = P7; P31 = P8 ; P21 =P9 ;
%     P = [P11 P12 P13 ;     P21 P22 P23 ;     P31 P32 P33] ;   
    
else
    F11 = F1;  F22 = F2 ;  F12 = F3 ; F21 = F4 ;
    F = [F11 F12 ;          F21 F22] ;
    
%     
%      P11 = P1;  P22 = P2 ;  P12 = P3 ; P21 = P4 ;
%     P = [P11 P12 ;          P21 P22] ;
end


detF = det(F) ; 



%%%%%


% Equivalence
nF = 9 ;
Fequiv = cell(nF,2) ;
for i=1:nF
    Fequiv{i,1} = ['F',num2str(i)];
    Fequiv{i,2} = ['FgradST(FROWS{',num2str(i),'})'];
end
% nstrain = 6 ;
% Cequiv = cell(nstrain,nstrain,2) ;
% for istrain =1:nstrain
%     for jstrain = 1:nstrain
%         Cequiv{istrain,jstrain,1} = ['C',num2str(istrain),num2str(jstrain)] ;
%         Cequiv{istrain,jstrain,2} = ['celastST(CROWS{',num2str(istrain),'},',num2str(jstrain),')'] ;
%     end
% end

diary('detFfile.txt')

%for i=1:size(CmatFIN,1)
 %   for j=i:size(CmatFIN,2)
  %      label_i = ['FROWS{',num2str(i),'}'] ;
   %     label_j = num2str(j) ;
    %    C_loc= ['celasLARGEmat(',label_i,',',label_j,') = '] ;
     %   Cstr = sym2str(CmatFIN(i,j)) ;
        
        %    if  j >=i
        detF = sym2str(detF) ; 
        for ireplace = 1:size(Fequiv,1)
            detF = strrep(detF,Fequiv{ireplace,1},Fequiv{ireplace,2}) ;
        end
%         for ireplace = 1:size(Cequiv,1)
%             for jreplace= 1:size(Cequiv,2)
%                 Cstr = strrep(Cstr,Cequiv{ireplace,jreplace,1},Cequiv{ireplace,jreplace,2}) ;
%             end
%         end
        
        %   else
        %                 Cstr= ['celasLARGEmat(',label_j,',',label_i,') = '] ;
        
        %  end
        disp(['detF =',detF,';']) ;
        
%     end
% end

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


