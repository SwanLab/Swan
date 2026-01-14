function V = SVDT_norm(V,Mchol,DATAIN)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Now we apply the weighted SVD over V so that they become
% mutually  M-orthogonal 
Xbar = Mchol*V ; 
DATAIN = DefaultField(DATAIN,'TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES',1e-6) ; 
TOL = DATAIN.TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES ;
if isempty(TOL)
    TOL = 0 ; 
end
DATALOC.RELATIVE_SVD = 1; 
[Ubar,S,Vbar] = SVDT(Xbar,TOL,DATALOC) ;
V = Mchol\Ubar ;

% Sjump  = zeros(size(S)); 
% Sjump(1) = 1; 
% Sjump(2:end) = S(2:end)./S(1:end-1) ;  
