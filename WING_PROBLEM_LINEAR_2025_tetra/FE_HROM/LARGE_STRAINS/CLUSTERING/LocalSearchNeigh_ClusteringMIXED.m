function [indCANDIDATES,sortedDISTANCES,indDIST ]=  ...
    LocalSearchNeigh_ClusteringMIXED(ipoints,d,nd,SNAPdisp,IndRest,TolStopSearch,MaxNumberNeighbor,DISTANCES_POINTS,...
    indGlobalSelf,MinNumberNeighbors_interTRAJ)

if nargin == 0
    load('tmp1.mat')
    DISTANCES_POINTS = [] ; 
end

% Compute the distance between this point and the remaining ones
if isempty(DISTANCES_POINTS)
DIST = bsxfun(@minus,SNAPdisp(:,IndRest),d) ;
nDIST  = sqrt(sum(DIST.^2,1)) ; % Euclidean distance
else
    nDIST = DISTANCES_POINTS(IndRest,ipoints)'; 
end


[sortedDISTANCES,indDISTloc] = sort(nDIST)  ;
indDIST =  IndRest(indDISTloc); %

% We seek to express d as a linear combination of its neighboring
% points
jpoints = 1;
ERROR_approx = 1e20 ;



while (jpoints <=MaxNumberNeighbor  && ERROR_approx > TolStopSearch)  ||   jpoints  <= MinNumberNeighbors_interTRAJ
    indCANDIDATES = [indGlobalSelf(:);indDIST(1:jpoints)] ;
    BASIS =  orth(SNAPdisp(:,indCANDIDATES))  ;
    ERROR_approx = norm(d - BASIS*(BASIS'*d))/nd ;
    jpoints = jpoints + 1;
end

indCANDIDATES = indCANDIDATES(length(indGlobalSelf)+1:end) ; 


% if   jpoints == MaxNumberNeighbor+1
%  %   disp(['Not converged ....Approximation error =',num2str(ERROR_approx)]) ;
% %    error(['Not converged for point =',num2str(ipoints)])
%     
% end

