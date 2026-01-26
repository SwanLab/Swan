function [VAR,celastST,FgradST,detFgrad,KcGEOunassemGLOloc,PdownsRBcoupROTc_i] = ...
    ResidualStaticCOROTfinal(OPERFE,D_QrotALL,VAR,LboolCallQ,MATPRO,VARint_n,DATA,dCqLOC,Qrot)
% Nodal residual, static, EIFEM, corotational
% See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 19-feb-2025, WEDNESDAY, 11:35, Balmes 185,  Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf

    
    % -----------------------------------------------------------------
    % Difference between total displacement (in local coordinates) and
    % rotational displacements 
    % ------------------------------------------------------------------
    % \dClocINCRE = \LboolCallQ    \dC   - \dCqLOC
    VAR.dClocINCRE = LboolCallQ*VAR.DISP -dCqLOC;      
    % --------------------------------------------------
    % UNASSEMBLED INTERNAL FORCES  (local coordinates)
    % --------------------------------------------------
    % LINEAR COMPONENT 
    %{(\FintCunassembLOC)^{lin}} =  \DiagC{\KcLINloc} \dClocINCRE
    FintCunassembLOC_lin = OPERFE.D_KcLINloc*VAR.dClocINCRE ; 
    
    % NONLINEAR COMPONENT  
    [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR_COROT_LRss(OPERFE,VAR,MATPRO,DATA,VARint_n) ; % OUTPUT VAR.PK1STRESS_incre
    %  \FintCunassembLOC)^{non}  = \DiagC{\BmatIst}^T \DiagC{\Wecm}  (\PoneFlocST -\PoneFlocLINst )
    FintCunassembLOC_non = OPERFE.D_BmatIst'*(OPERFE.D_Wecm)*VAR.PK1STRESS_incre ; 
    % TOTAL (linear + nonlinear)
    %  \FintCunassembLOC =  \ {(\FintCunassembLOC)^{lin}} +  \ {(\FintCunassembLOC)^{nonlin}}
    FintCunassembLOC = FintCunassembLOC_lin + FintCunassembLOC_non ;  
    % MODIFIED INTERNAL FORCES (TO ACCOUNT FOR THE COMPATIBILITY RESIDUAL)
    %  \KcGEOunassemGLOloc \defeq   \DiagC{\AspMAT}  \DiagC{\FintCunassemb} \DiagC{\PdownsRBcoupROTc}
    [KcGEOunassemGLOloc,PdownsRBcoupROTc_i] = KstiffGEO_unass_compute(OPERFE,DATA,D_QrotALL,FintCunassembLOC,Qrot) ; 
    % FintC_geo = LboolCall'*\KcGEOunassemGLOloc*dClocINCRE ;
     FintC_geo = OPERFE.LboolCall'*(KcGEOunassemGLOloc*VAR.dClocINCRE) ; 
    % ----------------------------
    % ASSEMBLED INTERNAL FORCES  (GLOBAL COORDINATES)
    % ----------------------------
    VAR.FINT = LboolCallQ'*FintCunassembLOC  +FintC_geo  ; 
    VAR.RESID  = VAR.FINT- VAR.FEXT;  
 