function R = OrthogonalComplement(M,Mchol,R,V,DATAIN)

PG = (V'*M*V)  ;
R = R  - V*(PG\(V'*M*R)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Now we apply the weighted SVD over R so that they become
% mutually  M-orthogonal 
Xbar = Mchol*R ; 
DATAIN = DefaultField(DATAIN,'TOL',1e-6) ; 
TOL = DATAIN.TOL ;
if isempty(TOL)
    TOL = 0 ; 
end
DATALOC.RELATIVE_SVD = 1; 
[Ubar,S,Vbar] = SVDT(Xbar,TOL,DATALOC) ;
R = Mchol\Ubar ;

% Sjump  = zeros(size(S)); 
% Sjump(1) = 1; 
% Sjump(2:end) = S(2:end)./S(1:end-1) ;  
