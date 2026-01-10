function [qRB,qDEF,rDEF] = DomainEquationAssemblyMULTI(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN,DATAwmethod,rRBglo)
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% EQUATIONS
%-------------------------
%------------------------
DATAIN = DefaultField(DATAIN,'SchurComplementSolution',0) ;
DATAIN = DefaultField(DATAIN,'UnconstrainedMinimization',0) ;


if DATAIN.SchurComplementSolution == 0  & DATAIN.MinimizationBoundaryWork == 0 & DATAIN.UnconstrainedMinimization==0
    
    
    if  DATAIN.ConstrainedWithReactionsAtInterfaces == 0
        DATAIN =DefaultField(DATAIN,'SimpleProjection',0);

        if DATAIN.SimpleProjection ==0 
        [qRB,qDEF,rDEF] = AssemblyStandard(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
            KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN);
        else
            [qRB,qDEF,rDEF] = AssemblySimpleProjection(nDOM,DATAwmethod,INDrig,INDdef,...
            KdomREDglo,Hqr,Fbar);
        end
        
    elseif  DATAIN.ConstrainedWithReactionsAtInterfaces == 1
        [qRB,qDEF,rDEF] =DomEqAsseConstrReaction(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
            KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN,DATAwmethod) ;
    else
        error('Option not implemented')
    end
    
    
    
elseif DATAIN.SchurComplementSolution == 1
    
    [qRB,qDEF,rDEF] =  DomainEquationAssemblySCHUR(nDOM,nRB,nDEF,nREAC,Pcomp, bCOMP,INDrig,INDdef,...
        KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN) ;
elseif  DATAIN.MinimizationBoundaryWork == 1
    
    [qRB,qDEF,rDEF] =  DomainEquationAssemblyENERGY1(nDOM,DATAwmethod,INDrig,INDdef,...
        KdomREDglo,Hqr,Fbar,DATAIN) ;
    
elseif  DATAIN.MinimizationBoundaryWork == 2
    
   [qRB,qDEF,rDEF] = AssemblyStandardMINBWORK(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
            KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAwmethod,rRBglo);
    
elseif  DATAIN.UnconstrainedMinimization == 1
    [qRB,qDEF,rDEF]  =    AssemblyUnconstrainedMinimization(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
        KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN) ;
    
else
    
    error('Option not implemented')
    
end