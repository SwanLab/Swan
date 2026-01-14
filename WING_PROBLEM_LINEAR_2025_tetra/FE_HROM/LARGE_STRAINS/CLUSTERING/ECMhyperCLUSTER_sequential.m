function [ECMdata_cluster,setCandidates,maX_INTERNAL_FORCES] = ECMhyperCLUSTER_sequential(BasisU_cluster,BasisPone_cluster,BasisStwo_cluster,...
    OPERFE,DISP_CONDITIONS,DATA,DATAoffline)

if nargin  == 0
    load('tmp2.mat')
end

diary off

nclusters = length(BasisPone_cluster) ;


disp('SEQUENTIAL METHOD FOR HYPERREDUCTION')
disp('********************************************************************')
disp(['Determining basis matrices for internal forces for each cluster'])
disp('********************************************************************')
BasisFint_cluster = cell(size(BasisStwo_cluster));
DOFl = DISP_CONDITIONS.DOFl ;
DOFr = DISP_CONDITIONS.DOFr ;

if isstruct(BasisU_cluster)
    Bst_DOFl = OPERFE.Bst(:,DOFl)*BasisU_cluster.BASIS ;
else
    Bst_DOFl = [] ;
end
maX_INTERNAL_FORCES= 0 ; 
for icluster = 1:nclusters
    if isstruct(BasisU_cluster)
        BstRED_l = Bst_DOFl*BasisU_cluster.coeff{icluster} ;
    else
        BstRED_l =OPERFE.Bst(:,DOFl)*BasisU_cluster{icluster} ;
    end
    %     if icluster == 553
    %         disp('borrar esto')
    %     end
    
    BasisFint_cluster{icluster} = QbasisMatrixIntegrand(BstRED_l,BasisPone_cluster{icluster},DATA,OPERFE.wSTs,DATAoffline) ;
    disp(['CLUSTER = ',num2str(icluster),';  N modes FINT =  ',num2str(size(BasisFint_cluster{icluster},2))]) ;
    maX_INTERNAL_FORCES = max(maX_INTERNAL_FORCES,size(BasisFint_cluster{icluster},2)) ;
end

% Sort the number of columns of each BasisFint_cluster
DATAoffline = DefaultField(DATAoffline,'OrderSequentialECM','SORTED') ; % If =2,
[dummy,nmodesFINT] = cellfun(@size,BasisFint_cluster) ;


DATAoffline = DefaultField(DATAoffline,'USE_GLOBAL_SUBSPACE_FOR_INTERNAL_FORCES',0) ;   % 6-Feb-2024
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/PAPER_SAW_ECM.mlx



if  DATAoffline.USE_GLOBAL_SUBSPACE_FOR_INTERNAL_FORCES  == 1
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/PAPER_SAW_ECM.mlx
    
    [ECMdata,~] = ECMcluster_global(BasisFint_cluster,...
        OPERFE,DATA,DATAoffline) ;
    setElements = large2smallREP(ECMdata.setPoints,DATA.MESH.ngaus) ;
    disp('****************************+')
    disp(['List of selected m = ',num2str(length(setElements)),' elements (for all clusters)'])
    disp(num2str(setElements'))
    setCandidates =ECMdata.setPoints;
    
    ECMdata_cluster =cell(size(BasisFint_cluster)) ;
    
    ECMdata_cluster(:) =  {ECMdata} ;
    
    
    setIndices = small2large(ECMdata.setPoints,DATA.MESH.nstrain) ;
    
    
    for icluster = 1:length(ECMdata_cluster)
        
        if isstruct(BasisStwo_cluster)
            BasisStwo= BasisStwo_cluster.BASIS*BasisStwo_cluster.coeff{icluster} ; ;   % Bssis matrix provided in a "multiscale" fashion
        else
            BasisStwo = BasisStwo_cluster{icluster} ;
        end
        BasisStwoZ = BasisStwo(setIndices,:) ;
        ECMdata_cluster{icluster}.coeff_reconstr_STWO =  (BasisStwoZ'*BasisStwoZ)\BasisStwoZ' ;
    end
    
    
    
else
    
    
    switch DATAoffline.OrderSequentialECM
        case 'Sequential'
            
            switch   DATAoffline.METHOD_PARTITION_SNAPSHOTS
                case  'SEQUENTIAL_unitrajectory_3'
                    NEW_ORDER_clusters = 1:length(nmodesFINT) ;
                otherwise
                    error('Option  not valid  for this type of partition')
            end
            
            
        case 'SORTED'
            [III,NEW_ORDER_clusters] = sort(nmodesFINT,'descend') ;
        case 'RANDOM'
            DATAoffline = DefaultField(DATAoffline,'NUMBER_RANDOM_RUNS_OrderSequentialECM',1)  ;
            
            TOTAL_number_perm =factorial(length(BasisFint_cluster));
            
            if  TOTAL_number_perm <= DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM
                NEW_ORDER_clusters = perms(1:length(BasisFint_cluster)) ;
            else
                NEW_ORDER_clusters = zeros(DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM,length(BasisFint_cluster)) ;
                for  irandom = 1:DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM
                    NEW_ORDER_clusters(irandom,:) = randperm(length(BasisFint_cluster)) ;
                end
                
            end
        case 'RANDN_fixed'
            rng(1)
            NEW_ORDER_clusters = randperm(length(BasisFint_cluster)) ;
            %  NEW_ORDER_clusters = 1:length(BasisFint_cluster) ;
        case 'PRINCIPAL_ANGLES_SUBSPACES'
            % Criterion based on the principla angles formed by the column
            % spaces of BasisFint_cluster
            NEW_ORDER_clusters  =  SortECM_sequentialANGLES(BasisFint_cluster) ;
    end
    
    % Next we run a loop over all clusters, starting with the cluster with
    % higher number of modes (i.e., points)
    disp('*************************************************++')
    disp(['Sequential ECM ...'])
    
    %ntrials = 2;
    setCandidates= [] ;
    ECMdata_cluster = [] ;
    lengthSETCAND = 1e+40 ;
    
    NumberCandidats= zeros(size(NEW_ORDER_clusters,1),1) ;
    
    itrial_BEST = 1;
    
    LIST_OF_CANDIDATES_BEST =[] ;
    for itrials = 1:size(NEW_ORDER_clusters,1)
        [ECMdata_cluster_test,setCandidates_test,LIST_OF_CANDIDATES_test] = ECM_local_clusteringLOOP(nclusters,NEW_ORDER_clusters(itrials,:),BasisStwo_cluster,OPERFE,DATAoffline,BasisFint_cluster,...
            DATA,BasisU_cluster,BasisPone_cluster)  ;
        if lengthSETCAND > length(setCandidates_test)
            setCandidates = setCandidates_test ;
            ECMdata_cluster = ECMdata_cluster_test ;
            lengthSETCAND = length(setCandidates) ;
            itrial_BEST = itrials ;
            LIST_OF_CANDIDATES_BEST =LIST_OF_CANDIDATES_test ;
        end
        NumberCandidats(itrials) =  length(setCandidates) ;
    end
    
    NEW_ORDER_clusters_best = NEW_ORDER_clusters(itrial_BEST,:) ;
    save('NEW_ORDER_clusters_best_ALL.mat','NumberCandidats','NEW_ORDER_clusters_best','ECMdata_cluster','setCandidates','LIST_OF_CANDIDATES_BEST') ;
    
    disp('***********************************************************************+')
    disp(['Number of ECM points per random run = ',num2str(NumberCandidats')])
    disp('***********************************************************************+')
    save('NumberCandidatsRANDOM.mat','NumberCandidats')
    disp(['Total number of points (had they been  selected independently) = ',num2str(sum(nmodesFINT))])
    disp(['Total number of ECM candidate points = ',num2str(length(setCandidates))])
    setElements = large2smallREP(setCandidates,DATA.MESH.ngaus) ;
    disp('****************************+')
    disp(['List of selected m = ',num2str(length(setElements)),' elements (for all clusters)'])
    disp(num2str(setElements'))
    %   clipboard('copy',num2str(setElements'));
    
end