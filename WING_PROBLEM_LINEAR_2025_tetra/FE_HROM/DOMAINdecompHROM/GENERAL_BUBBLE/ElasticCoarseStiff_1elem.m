function KcELAS = ElasticCoarseStiff_1elem(EIFEoper_all,indCHOSEN,TRANSF_COORD,ndofsBUB)
if nargin == 0
    load('tmp1.mat')
end
        % STIFFNESS MATRIX COARSE, elastic
        %---------------------------------
        % 5-March-2024
        % See explanation in 
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        % JAHO, Balmes 185, Barcelona
      
% % Block boundary-boundary DOFs
DOFsB =  EIFEoper_all(indCHOSEN).INFO.DOFsB ;  % DOFs boundary
DOFsBUB =  EIFEoper_all(indCHOSEN).INFO.DOFsBUB ;  % DOFs bubble
% Block matrices 
KcBBelL= EIFEoper_all(indCHOSEN).Kcoarse(DOFsB,DOFsB) ;
KcBBUBelL= EIFEoper_all(indCHOSEN).Kcoarse(DOFsB,DOFsBUB) ;
KcBUBBelL = EIFEoper_all(indCHOSEN).Kcoarse(DOFsBUB,DOFsB) ;
KcBUBBUBelL= EIFEoper_all(indCHOSEN).Kcoarse(DOFsBUB,DOFsBUB) ;
 
% INTRODUCE ROTATIONS (might not be very efficient, check it)
% LEFT
nDOFSloc =  (size(KcBBelL,1)) ;  % 
nSD = length(TRANSF_COORD.TRANSLATION); % Number of spatial dimensions
nnodesLOC = nDOFSloc/nSD ;  % Number of integration points
ROTglo = cell(1,nnodesLOC) ;
ROTglo(:) = {sparse(TRANSF_COORD.ROTATION)} ;
ROTglo =  blkdiag(ROTglo{:}) ; 
ROTgloT = cell(1,nnodesLOC) ;
ROTgloT(:) = {sparse(TRANSF_COORD.ROTATION')} ;
ROTgloT =  blkdiag(ROTgloT{:}) ;   
% 
% \ {\KcBBel}  = (\lambda^{\nSD-2})  \Qall_n  \KcBBelL  \Qall_n^T
lambda_n2 = TRANSF_COORD.SCALEFACTOR^(nSD-2) ; 
KcBBel = lambda_n2*ROTglo*KcBBelL*ROTgloT; 
%  \ {\KcBBUBel} = (\lambda^{\nSD-2})  \Qall_n  \KcBBUBelL  
KcBBUBel = lambda_n2*ROTglo*KcBBUBelL;
%  \ {\KcBUBBel} = (\lambda^{\nSD-2})   \KcBUBBelL   \Qall_n^T
KcBUBBel = lambda_n2*KcBUBBelL*ROTgloT; 
%  \{\KcBUBBUBel} = (\lambda^{\nSD-2})   \KcBUBBUBelL    
KcBUBBUBel = lambda_n2*KcBUBBUBelL; 


KcELAS = ExtendedKmatrixBUB(KcBBel,KcBBUBel,KcBUBBel,KcBUBBUBel,nSD,ndofsBUB) ; 

 