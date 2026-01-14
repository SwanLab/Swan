function [RESIDUALF,DATAINM ]= ComputeResidual(qDEF,rDEF,rRB,fextDOMglo,BasisUdef,BasisRdef,BasisRrb,K,DATAINM);

if nargin == 0
    load('tmp.mat')
end

DATAINM = DefaultField(DATAINM,'RESIDUAL_FORCES',[]) ; 
DATAINM.RESIDUAL_FORCES = DefaultField(DATAINM.RESIDUAL_FORCES,'COMPUTE',1) ; 

RESIDUALF.ALL = [] ; 
RESIDUALF.DOMAIN_AVGNODE = [] ; 

nDEF = size(BasisUdef{1},2) ; 
nDEFr = size(BasisRdef{1},2) ; 
nRB = size(BasisRrb,2) ;
ndom = length(fextDOMglo) ;
Kmax = max(max(K)) ; 
K = K*BasisUdef{1} ; 
if nRB ==6
    ndim = 3; 
else
    ndim = 2; 
end

DATAINM  = DefaultField(DATAINM,'FACTOR_MULTIPLY_RESIDUAL_FORCES',1) ;  
 

if DATAINM.RESIDUAL_FORCES.COMPUTE ==1 
    RESIDUALF.ALL = cell(ndom,1) ; 
    RESIDUALF.NORM_AVERAGE = zeros(ndom,1) ; 
        RESIDUALF.MAX_NORM = zeros(ndom,1) ; 

    disp('Computing residual forces')
    ifinU = 0 ; ifinRdef = 0 ; ifinRrb = 0 ; 
    for idom = 1:ndom
        iniU = ifinU + 1 ; ifinU = ifinU + nDEF ; 
        iniRdef = ifinRdef + 1 ; ifinRdef = ifinRdef + nDEFr ;
        iniRrb = ifinRrb + 1 ; ifinRrb = ifinRrb + nRB ;
        RES = K*qDEF(iniU:ifinU) - fextDOMglo{idom} + (BasisRrb*rRB(iniRrb:ifinRrb)+ BasisRdef{1}*rDEF(iniRdef:ifinRdef)   ) ;
        RESIDUALF.ALL{idom} = RES; 
        nRES_pernode= norm(RES)/(length(RES)/ndim) ;
        RESIDUALF.NORM_AVERAGE(idom) = nRES_pernode ; 
        RES = reshape(RES,ndim,[]) ; 
        RES = sqrt(sum(RES.^2,1));
        RESIDUALF.MAX_NORM(idom) = max(RES) ;
        
    end
    RESIDUALF.ALL = cell2mat(RESIDUALF.ALL) ; 
    disp('...Done')
    
    RESIDUALF.ALL = RESIDUALF.ALL*DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES ;
    RESIDUALF.MAX_NORM = RESIDUALF.MAX_NORM*DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES ;
    RESIDUALF.NORM_AVERAGE = RESIDUALF.NORM_AVERAGE*DATAINM.FACTOR_MULTIPLY_RESIDUAL_FORCES ;

end





% figure(178)
% hold on 
% xlabel('Domain')
% ylabel('Max.Norm.ResidF')
% plot(RESIDUALF.MAX_NORM)
% figure(179)
% hold on 
% xlabel('Domain')
% ylabel('Resid. Per Node')
%  plot(RESIDUALF.NORM_AVERAGE,'r')




 
