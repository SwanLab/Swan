function [Fbar,fextDOMglo,reactDOMrbGLO] = ...
    ExternalForcesReduced(reactDOMrbGLO,fextDOMglo,BasisUdefGLO,DATAINM,BasisUrbGLO,INDrig,INDdef)

if nargin == 0
    load('tmp.mat')
end
reactDOMrbGLO = cell2mat(reactDOMrbGLO) ;
reactDOMrbGLO = (reactDOMrbGLO(:)) ;
fextDOMglo  = cell2mat(fextDOMglo) ;
fextDOMglo = (fextDOMglo(:)) ;
if DATAINM.MinimizationBoundaryWork == 0 | DATAINM.MinimizationBoundaryWork == 2
    
    
    Fbar = BasisUdefGLO'*(fextDOMglo+reactDOMrbGLO) ;
    
elseif DATAINM.MinimizationBoundaryWork == 1
    nmodes =length(INDrig) + length(INDdef) ;
    Fbar = zeros(nmodes,1) ;
    Fbar(INDdef) = BasisUdefGLO'*fextDOMglo ;
    Fbar(INDrig) = BasisUrbGLO'*fextDOMglo    ;
    
else
    error('Option not implemented')
    
end