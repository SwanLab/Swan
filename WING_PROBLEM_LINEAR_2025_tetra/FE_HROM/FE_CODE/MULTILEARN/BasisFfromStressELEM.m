function [BasisF,SingVal_intf,ngausLOC,nstrain] = ...
    BasisFfromStressELEM(BasisS,SingVal_stress,BdomRED, Wdom,DATAIN,nelem)
% Similar to BasisFfromStress.m, but summing up the contribution of each
% element 
%dbstop('4')
if nargin ==0
    load('tmp.mat')
   % DATACUB = [] ;
end
DATAIN = DefaultField(DATAIN,'DATACUB',[]) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ; 
DATACUB = DATAIN.CUBATURE ; 
DATACUB= DefaultField(DATACUB,'IncludeSingularValuesStress',1) ; % Include Singular values stresses
if DATACUB.IncludeSingularValuesStress == 1
    BasisS = bsxfun(@times,BasisS',SingVal_stress)' ;
end
ngaus = length(Wdom);
nstrain = size(BdomRED,1)/length(Wdom); 
nDEF = size(BdomRED,2) ; % Number of displacement modes
nBASISs = size(BasisS,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
for I = 1:nDEF
    SNAPloc = zeros(ngaus,nBASISs) ;
    for istrain = 1:nstrain
        B_stress = bsxfun(@times,BasisS(istrain:nstrain:end,:),BdomRED(istrain:nstrain:end,I)); 
     %   SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ; 
       SNAPloc = SNAPloc + bsxfun(@times,B_stress,Wdom) ; 
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;     
end


% Summing up the contribution of each element 

SNAPforceS_elem = zeros(nelem,size(SNAPforceS,2));
ngausLOC  = ngaus/nelem ; 

for igaus = 1:ngausLOC
    SNAPforceS_elem = SNAPforceS_elem + SNAPforceS(igaus:ngausLOC:end,:) ; 
end


SNAPforceS = SNAPforceS_elem  ; 



if ~exist('RSVDT')
    addpath('SVDlibrary')
end




hold on
nfigure = 100;
LEGENDG = 'Internal Virtual Work ' ;
COLOR = 'k-' ;
%dbstop('40')
DATAIN = DefaultField(DATAIN,'TOL_LOC_InternalForces',1e-6) ; 
DATAIN.TOL_LOC = DATAIN.TOL_LOC_InternalForces  ; 
[U,S,V,h3,h4] = SVD_and_error(SNAPforceS,nfigure,LEGENDG,[],COLOR,DATAIN ) ;
hold on
legend([h3 ],{'Int. virt. work (element assembly)'})
legend([h4 ],{'Int. virt. work  (element assembly)'})

BasisF = U  ; % = [2] ;
SingVal_intf = S ;  
% save(DATA.NAMEWS,'-append','BasisS','SingVal_intf');


