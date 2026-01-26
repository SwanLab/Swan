function [DATAoffline,IND_loc,SNAPSHOTS_PER_CLUSTER] = OverlappingClusters(DATAoffline,DISTANCE_CENT,IND_loc,INDEX_SELECT)

if nargin == 0
    load('tmp.mat')
end


% For each snapshot, determine which is the two closest centroids

[aaa,CLOSEST_CENTROIDS] = sort(DISTANCE_CENT,2) ;

% Two closest centroids to each snapshot
DATAoffline = DefaultField(DATAoffline,'nneig_search',2) ; 
nneig_search  = DATAoffline.nneig_search;  

CLOSEST_CENTROIDS = CLOSEST_CENTROIDS(:,1:nneig_search) ;


nclusters = size(DISTANCE_CENT,2) ;

CNclusters = zeros(nclusters,nclusters) ;
ncluster_neig = zeros(nclusters,1) ;

% Loop over all snapshots.
for isnap = 1:size(CLOSEST_CENTROIDS,1)
    
    if  nneig_search == 2
        k = CLOSEST_CENTROIDS(isnap,1) ;
        l = CLOSEST_CENTROIDS(isnap,2) ;
        
        if ~any(CNclusters(k,1:ncluster_neig(k)) == l)
            ncluster_neig(k)  = ncluster_neig(k) + 1;
            CNclusters(k,ncluster_neig(k)) = l ;
        end
        
        if ~any(CNclusters(l,1:ncluster_neig(l)) == k)
            ncluster_neig(l)  = ncluster_neig(l) + 1;
            CNclusters(l,ncluster_neig(l)) = k ;
        end
        
        
    else 
        % General case 
        
        for  ineig = 1:nneig_search   % Loop over neighboring clusters
            k = CLOSEST_CENTROIDS(isnap,ineig) ;
            for  jneig = 1:nneig_search
                if  ineig ~= jneig
                    l = CLOSEST_CENTROIDS(isnap,jneig) ;                    
                    if ~any(CNclusters(k,1:ncluster_neig(k)) == l)
                        ncluster_neig(k)  = ncluster_neig(k) + 1;
                        CNclusters(k,ncluster_neig(k)) = l ;
                    end                    
                    if ~any(CNclusters(l,1:ncluster_neig(l)) == k)
                        ncluster_neig(l)  = ncluster_neig(l) + 1;
                        CNclusters(l,ncluster_neig(l)) = k ;
                    end                    
                end                
            end      
        end
        
        
        
        
    end
    
end

aaa = max(ncluster_neig) ;

CNclusters = CNclusters(:,1:aaa) ;




NameColums = cell(1,nclusters) ;
for i=1:length(NameColums)
    NameColums{i} = ['i=',num2str(i)] ;
end

NameRows = cell(1,aaa) ;
for i=1:length(NameRows)
    NameRows{i} = ['k=',num2str(i)] ;
end




sz = size(CNclusters);
varTypes = cell(size(NameColums)) ;
varTypes(:) = {'int'} ;
varNames = NameColums ;
T2 = array2table(CNclusters,'VariableNames',NameRows,'RowNames',varNames);
table2latex(T2,'Inter_clusterCONN.txt') ;
disp(T2)

DATAoffline.CNclusters_distance = CNclusters ;

DATAoffline = DefaultField(DATAoffline,'UseOverlapping',0);

SNAPSHOTS_PER_CLUSTER = cell(1,length(ncluster_neig)) ;

% From washabaugh2012nonlinear.pdf
for j = 1:length(ncluster_neig)
    Wj = find(IND_loc==j) ; %  Number of states in cluster  j
    nsnap = length(Wj);
    nneig =  ncluster_neig(j) ;     % number of its neighboring clusters
    nadd  = ceil(DATAoffline.UseOverlapping*nsnap/nneig) ;
    
    if DATAoffline.UseOverlapping == 0
        SNAPSHOTS_PER_CLUSTER{j} = INDEX_SELECT(Wj) ;
    else
        for iloc = 1:nneig
            i = CNclusters(j,iloc) ;  % Neighboring cluster
            % Find the Nadd,j states nearest to their common boundary and
            % add these to ÃWj
            indSNAP_i = find(IND_loc==i) ; % Indexes snapshots cluster i
            DIST_j = DISTANCE_CENT(indSNAP_i,j) ; % Distance states cluster i to centroid cluster j
            [DIST_j_sorted,IndexesSorted] = sort(DIST_j) ;
            nadd_i = min(length(IndexesSorted),nadd) ; % Number of neighboring clusters to add
            IndexesSortedLocal = IndexesSorted(1:nadd_i) ;
            IndexGlobal = indSNAP_i(IndexesSortedLocal) ; 
            Wj = [Wj;IndexGlobal]  ; 
        end
        
        SNAPSHOTS_PER_CLUSTER{j} = INDEX_SELECT(unique(Wj)) ; 
        
    end
    
end





%%

