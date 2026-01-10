clc
clear all
%--------------------------------------------------------------------------
% Script: SymbolicCauchyStressGeneration.m
%
% PURPOSE:
%   Symbolically compute the Cauchy stress tensor from the first Piola-Kirchhoff
%   stress tensor via:
%
%       σ = (1/det(F)) * Fᵗ * P
%
%   and export the result in MATLAB-readable syntax adapted for vectorized
%   implementations using pre-stored indices (e.g., FROWS, CROWS).
%
% FUNCTIONALITY:
%   - Defines symbolic deformation gradient F and first Piola-Kirchhoff tensor P.
%   - Computes σ = Fᵗ P / det(F)
%   - Extracts σ in Voigt form.
%   - Rewrites each component as a vectorized MATLAB assignment using
%     field references: FgradST(FROWS{i}) and PoneST(FROWS{i}).
%   - Writes the results to 'CuachyStress.txt' file for direct copy-paste
%     into `CauchyStressFromPK1.m`.
%
% DIMENSIONALITY:
%   - Works for both 2D and 3D (set via `ndim`).
%
% OUTPUT:
%   A series of MATLAB vectorized assignments of the form:
%
%     CauchyStress(CROWS{1}) = (FgradST(FROWS{1}).*PoneST(FROWS{1}) + ...) / detF;
%
% USAGE:
%   - Run once to generate or verify symbolic expressions.
%   - Modify if tensor-product structure changes or for different stress measures.
%
% LOCATION OF SYMBOLIC VERSION:
%   - See symbolic derivation also in:
%     /SYMBOLIC/SymCauchyStress.m
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 10-Dec-2020
%--------------------------------------------------------------------------


syms  F1 F2 F3 F4 F5 F6 F7 F8 F9
syms  P1 P2 P3 P4 P5 P6 P7 P8 P9 detF
ndim = 3;
if ndim == 3
    nF = 9 ; 
    nstrain = 6 ; 
    F11 = F1; F22 = F2 ; F33 = F3;  F23 = F4; F13 = F5; F12 = F6 ; F32 = F7; F31 = F8 ; F21 =F9 ;
    F = [F11 F12 F13 ;     F21 F22 F23 ;     F31 F32 F33] ;
    
    
    P11 = P1; P22 = P2 ; P33 = P3;  P23 = P4; P13 = P5; P12 = P6 ; P32 = P7; P31 = P8 ; P21 =P9 ;
    P = [P11 P12 P13 ;     P21 P22 P23 ;     P31 P32 P33] ;
    
    
    sigmaS = F.'*P/detF ;
    
    sigma{1} = sigmaS(1,1) ;
    sigma{2} = sigmaS(2,2) ;
    sigma{3} = sigmaS(3,3) ;
    sigma{4} = sigmaS(2,3) ;
    sigma{5} = sigmaS(1,3) ;
    sigma{6} = sigmaS(1,2) ;
    
    
    
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


