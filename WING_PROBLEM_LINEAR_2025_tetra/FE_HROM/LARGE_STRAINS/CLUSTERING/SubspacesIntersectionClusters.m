function [TransMatrix,CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,UPPERBOUND_matrix] = ...
    SubspacesIntersectionClusters(BasisU_cluster,DATAoffline,SNAPdisp,BasisU_ALL,DistanceCentroid)

if nargin == 0 
    load('tmp.mat')
end

NameColums = cell(1,length(BasisU_cluster)) ;
for i=1:length(NameColums)
    NameColums{i} = ['i=',num2str(i)] ;
end

 
NEW_METHOD = 1 ; % More efficient, in terms of CPU speed,  23th July-2022 

if NEW_METHOD ==1

    if DATAoffline.EfficentMemoryOption == 0
        [CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,TransMatrix,UPPERBOUND_matrix,DIMENSION_INTERSECTION_SPACEneigh] ...
            = SubspacesIntersectionLCnew(BasisU_cluster,DATAoffline,DistanceCentroid) ;
    else
        % This new routine alleviates the computation burden of
        % constructing NCLUSTER x NCLUSTER cells (23th July 2022)
         [CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,TransMatrix,UPPERBOUND_matrix,DIMENSION_INTERSECTION_SPACEneigh] ...
            = SubspacesIntersectionLCeff(BasisU_cluster,DATAoffline,DistanceCentroid) ;
    end
else
    [CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,TransMatrix,UPPERBOUND_matrix,DIMENSION_INTERSECTION_SPACEneigh] ...
    = SubspacesIntersectionLOCclust(BasisU_cluster,DATAoffline,DistanceCentroid) ; 
end



DATAoffline = DefaultField(DATAoffline,'CNclusters_distance',[]) ; 

if ~isempty(DATAoffline.CNclusters_distance)
% LOWERBOUND_matrix =  (LOWERBOUND_matrix+LOWERBOUND_matrix') ;

disp('Lower bound *********************+')

format short

sz = size(LOWERBOUND_matrix);
varTypes = cell(size(NameColums)) ;
varTypes(:) = {'double'} ;
varNames = NameColums ;
T2 = array2table(LOWERBOUND_matrix,'VariableNames',varNames,'RowNames',varNames);
table2latex(T2,'LowerBoundErrorTrans.txt') ;
disp(T2)


%%%%% 
if ~isempty(DATAoffline.CNclusters_distance)
TRANS_ERROR_NEIGH = zeros(size(DATAoffline.CNclusters_distance)) ; 
end

for iclusters = 1:size(DATAoffline.CNclusters_distance,1) 
   iloc = find(DATAoffline.CNclusters_distance(iclusters,:) > 0) ; 
   TRANS_ERROR_NEIGH(iclusters,iloc)  = LOWERBOUND_matrix(iclusters,iloc) ; 
end

NameRows = NameColums; 



NameColumns = cell(1,size(TRANS_ERROR_NEIGH,2)) ;
for i=1:length(NameColumns)
    NameColumns{i} = ['n=',num2str(i)] ;
end
disp(['Lower bound int. space'])

sz = size(TRANS_ERROR_NEIGH);
varTypes = cell(size(NameColums)) ;
varTypes(:) = {'double'} ;
varNames = NameColumns ;
T2 = array2table(TRANS_ERROR_NEIGH,'VariableNames',varNames,'RowNames',NameRows);
table2latex(T2,'TransErrorNeigh.txt') ;
disp(T2)

disp(['Dimension intersection space neigh. for TOL=',num2str(DATAoffline.LOWERBOUND_error_transition_neigh)])
varTypes(:) = {'int'} ;
varNames = NameColumns ;
T2 = array2table(DIMENSION_INTERSECTION_SPACEneigh,'VariableNames',varNames,'RowNames',NameRows);
table2latex(T2,'DimIntSpace.txt') ;
disp(T2)



end





% format long g 
% 
% disp('Upper bound *********************+')
% 
% format short
% 
% sz = size(UPPERBOUND_matrix);
% varTypes = cell(size(NameColums)) ;
% varTypes(:) = {'double'} ;
% varNames = NameColums ;
% T2 = array2table(UPPERBOUND_matrix,'VariableNames',varNames,'RowNames',varNames);
% table2latex(T2) ;
% disp(T2)
% 
% format long g 

% 
 