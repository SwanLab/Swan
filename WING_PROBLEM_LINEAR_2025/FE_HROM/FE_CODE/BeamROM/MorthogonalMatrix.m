function [BasisREFERENCE,S,IS_PARALLEL] = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface,DATAIN)

nb = norm(BasisREFERENCE,'fro') ; 

PG = (BeamModesInterface'*M*BeamModesInterface)  ;
BasisREFERENCE = BasisREFERENCE  - BeamModesInterface*(PG\(BeamModesInterface'*M*BasisREFERENCE)) ;

nb_after = norm(BasisREFERENCE,'fro') ; 

mag_projection = nb_after/nb ; 

if mag_projection <1e-10 
    warning('No orthogonal complement here: the basis is contained in the rigid body space')
    IS_PARALLEL = 1; 
else
    IS_PARALLEL = 0 ; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Now we apply the weighted SVD over BasisREFERENCE so that they become
% mutually  M-orthogonal 
Xbar = Mchol*BasisREFERENCE ; 
DATAIN = DefaultField(DATAIN,'TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES',1e-6) ; 
TOL = DATAIN.TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES ;
if isempty(TOL)
    TOL = 0 ; 
end
DATALOC.RELATIVE_SVD = 1; 
[Ubar,S,Vbar] = SVDT(Xbar,TOL,DATALOC) ;
BasisREFERENCE = Mchol\Ubar ;


 

% Sjump  = zeros(size(S)); 
% Sjump(1) = 1; 
% Sjump(2:end) = S(2:end)./S(1:end-1) ;  
