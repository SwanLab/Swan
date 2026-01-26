function [CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,...
    TransMatrix,UPPERBOUND_matrix,DIMENSION_INTERSECTION_SPACEneigh] = ...
    SubspacesIntersectionLCeff(BasisU_cluster,DATAoffline,DistanceCentroid)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/09_HOMOG_2D_1T/MakingItEfficient.mlx
if nargin == 0
    load('tmp.mat')
end

%% ANALYZING INTERSECTION BETWEEN SUBSPACES
% *****************************************
LOWERBOUND_matrix = sparse(length(BasisU_cluster),length(BasisU_cluster)) ;
UPPERBOUND_matrix = sparse(length(BasisU_cluster),length(BasisU_cluster)) ;

DATAoffline = DefaultField(DATAoffline,'LOWERBOUND_error_transition',1e-5); % Lower bound
DATAoffline = DefaultField(DATAoffline,'LOWERBOUND_error_transition_neigh',1e-7); % Lower bound


 %  TransMatrix
% clusters
% This matrix is only defined for neighboring clusters included in
% DistanceCentroid.DistanceCentroid.NEIGHBORING_CLUSTERS,
% The number of neighbors per cluster is given in DistanceCentroid.NumberOfNeighbors

% Accordingly, we make 
nmaxNEIGH = max(DistanceCentroid.NumberOfNeighbors) ; 
% and 
TransMatrix = cell(length(BasisU_cluster),nmaxNEIGH) ;




CONNECTIVITY_CLUSTERS =cell(length(BasisU_cluster),1);

DATAoffline = DefaultField(DATAoffline,'CNclusters_distance',[]) ;
if ~isempty(DATAoffline.CNclusters_distance)
    DIMENSION_INTERSECTION_SPACEneigh = zeros(size(DATAoffline.CNclusters_distance)) ;
else
    DIMENSION_INTERSECTION_SPACEneigh = [] ;
end

for i = 1:length(BasisU_cluster)
    NEIGH =  DistanceCentroid.NEIGHBORING_CLUSTERS{i}   ;
    for jloc = 1:length(NEIGH)  %
        j = NEIGH(jloc) ;
        %         if i ~= j
        %             ISneigh = 1; % It is a neighboring cluster
        %             if ~isempty(DistanceCentroid.NEIGHBORING_CLUSTERS)
        %                ISneigh =  ismember(j,DistanceCentroid.NEIGHBORING_CLUSTERS{i}) ;
        %             end
        %
        %             if ISneigh == 1
        %   disp(['i=',num2str(i),'; j=',num2str(j)]) ;
      %  TransMatrix{i,j} = BasisU_cluster{i}'*BasisU_cluster{j} ;
        
        TransMatrix{i,jloc} = BasisU_cluster{i}'*BasisU_cluster{j} ;
        [UU,SS,VV] = SVDT(TransMatrix{i,jloc}) ;
        LOWERBOUND =  real(sqrt(1-SS(1)^2)) ;
        LOWERBOUND_all =  real(sqrt(1-SS.^2)) ;
        
        [IndZEROall] = find(LOWERBOUND_all <= DATAoffline.LOWERBOUND_error_transition_neigh) ;
        
        %
        if ~isempty( DIMENSION_INTERSECTION_SPACEneigh )
            IJ =  find(j ==  DATAoffline.CNclusters_distance(i,:)) ;
            if ~isempty(IJ)
                DIMENSION_INTERSECTION_SPACEneigh(i,IJ) = length(IndZEROall) ;
            end
            
        end
        
        if size(TransMatrix{i,jloc},2) > size(TransMatrix{i,jloc},1)
            Upperbound = 1 ;
        else
            Upperbound = real(sqrt(1-SS(end)^2)) ;
        end
        
        %     disp(['LOWER BOUND ERROR = ',num2str(LOWERBOUND)]) ;
        [IndZERO] = find(LOWERBOUND <= DATAoffline.LOWERBOUND_error_transition) ;
        if ~isempty(IndZERO)
            CONNECTIVITY_CLUSTERS{i} = [CONNECTIVITY_CLUSTERS{i},j] ;
            % IntersectionMatrices{i,j} = UU(:,IndZERO) ;
        end
        LOWERBOUND_matrix(i,j)  = LOWERBOUND ;
        UPPERBOUND_matrix(i,j)  = Upperbound ;
        %      end
    end
end
 

