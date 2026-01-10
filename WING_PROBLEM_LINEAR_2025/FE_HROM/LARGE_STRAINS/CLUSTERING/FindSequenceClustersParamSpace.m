function CLUSTERS = FindSequenceClustersParamSpace(PATH_PARAM,DT,DATAoffline) ;
if nargin == 0
    load('tmp2.mat')
end
nsteps = size(PATH_PARAM,1) ; 
CLUSTERS =zeros(nsteps,1)  ; 
for istep = 1:nsteps
    POINT_param = PATH_PARAM(istep,:) ; % Parameter for time step istep
    POINT_refdomain = zeros(size(POINT_param)) ; 
    % Convert into reference domain
    for idim = 1:length(POINT_param)
        POINT_refdomain(:,idim) = (POINT_param(:,idim)-DT.BOUNDARIES_MIN(idim))/(DT.BOUNDARIES_MAX(idim)-DT.BOUNDARIES_MIN(idim)) ;
    end    
    % Find nearest point
    if size(DT.CN,2) >2
        CLUSTERS(istep) = tsearchn(DT.COOR,DT.CN,POINT_refdomain) ;
    else
        % 1D problem 
        DIST_cent = zeros(size(DT.CENTROID,1),1) ; 
        for idim = 1:length(POINT_param)
           DIST_cent  = DIST_cent  + (DT.CENTROID(:,idim)-POINT_refdomain(idim)).^2 ; 
        end
        
        [minDIST,IndicesCandidate] = sort(DIST_cent) ;
        
        if abs(minDIST(1)-minDIST(2)) <1e-16
            IndexNew = setdiff(IndicesCandidate(1:2),CLUSTERS(istep-1)) ; 
        else
            IndexNew = IndicesCandidate(1) ; 
        end
        CLUSTERS(istep) = IndexNew ; 
        
        
        
%         DT.CENTROID
%         CLUSTERS(istep) = knnsearch(DT.CENTROID,POINT_refdomain) ; 
    end    
end