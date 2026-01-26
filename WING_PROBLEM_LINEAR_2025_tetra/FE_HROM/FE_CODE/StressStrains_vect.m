function [strainGLOv,stressGLOv,DATAOUT,strainGLO,stressGLO,ngausLOC,wST] = ...
    StressStrains_vect(d,Bst,Cglo,DATA,wST,nelem,nstrain,ngaus,ndim)


strainGLOv = Bst*d ;
stressGLOv = (Cglo*strainGLOv)./wST ; % Change 6-Nov-2015

DATAOUT.strain = strainGLOv ;
DATAOUT.stress = stressGLOv ;


% For GID, we have to change the order of the components within the
% vector

if DATA.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 1
    [stressGLOv wSTelem] = AverageStressOnElements(stressGLOv,wST,nelem,nstrain,ngaus) ;  
    
    [strainGLOv] = AverageStressOnElements(strainGLOv,wST,nelem,nstrain,ngaus) ;
    ngausLOC = 1 ;
    wST = wSTelem ;
    [ DATAOUT.stressVONMISES ] =  VonMises_Stress(reshape(stressGLOv,nstrain,[] ))' ;
    
else
    ngausLOC = ngaus ;
end

ncomp = ngausLOC*nstrain*nelem  ;
strainGLO = strainGLOv;
stressGLO = stressGLOv;

if ndim == 3
    indGID = [1 2 3 6 4 5] ;
else
    if DATA.StrainStressWith4Components == 1
        indGID = [1 2 3 4] ;  % See GID's help
    else
        indGID = [1 2 3 ] ;
    end
end

for iglo = 1:nstrain
    strainGLO(iglo:nstrain:ncomp)  =  strainGLOv(indGID(iglo):nstrain:ncomp) ;
    stressGLO(iglo:nstrain:ncomp)  =  stressGLOv(indGID(iglo):nstrain:ncomp) ;
end


strainGLO =  reshape(strainGLO,ngausLOC*nstrain,[]) ;
stressGLO =  reshape(stressGLO,ngausLOC*nstrain,[]) ;

