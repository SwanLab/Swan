function [ECMdata_cluster,setCandidates] = LocalECM(wFE,A,DATAoffline)

 if nargin == 0
     load('tmp.mat')
 end

% Sort the number of columns of each A
DATAoffline = DefaultField(DATAoffline,'OrderSequentialECM','SORTED') ; % If =2,
[dummy,p_i] = cellfun(@size,A') ;

switch DATAoffline.OrderSequentialECM
    case 'Sequential'         
        switch   DATAoffline.METHOD_PARTITION_SNAPSHOTS
            case  'SEQUENTIAL_unitrajectory_3' 
                  NEW_ORDER_clusters = 1:length(p_i) ; 
            otherwise
                error('Option  not valid  for this type of partition') 
        end        
    case 'SORTED'
        [III,NEW_ORDER_clusters] = sort(p_i,'descend') ;
    case 'RANDOM'
        DATAoffline = DefaultField(DATAoffline,'NUMBER_RANDOM_RUNS_OrderSequentialECM',1)  ;
        NEW_ORDER_clusters = zeros(DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM,length(A)) ;
        for  irandom = 1:DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM
            NEW_ORDER_clusters(irandom,:) = randperm(length(A)) ;
        end
    case 'RANDN_fixed'
        rng(1)
        NEW_ORDER_clusters = randperm(length(A)) ;
        %  NEW_ORDER_clusters = 1:length(A) ;
    case 'PRINCIPAL_ANGLES_SUBSPACES'
        % Criterion based on the principla angles formed by the column
        % spaces of A
       NEW_ORDER_clusters  =  SortECM_sequentialANGLES(A) ; 
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

for itrials = 1:size(NEW_ORDER_clusters,1)
    
    [ECMdata_cluster_test,setCandidates_test] = LocalECMone(NEW_ORDER_clusters(itrials,:),wFE,DATAoffline,A)  ;
    
    if lengthSETCAND > length(setCandidates_test)
        setCandidates = setCandidates_test ;
        ECMdata_cluster = ECMdata_cluster_test ;
        lengthSETCAND = length(setCandidates) ;
    end
    
    NumberCandidats(itrials) =  length(setCandidates) ;
    
end

disp('***********************************************************************+')
disp(['Number of ECM points per random run = ',num2str(NumberCandidats')])
disp('***********************************************************************+')
save('NumberCandidatsRANDOM.mat','NumberCandidats')


disp(['Total number of points (had they been  selected independently) = ',num2str(sum(p_i))])

disp(['Total number of ECM candidate points = ',num2str(length(setCandidates))])

% 
% setElements = large2smallREP(setCandidates,DATA.MESH.ngaus) ;
% disp('****************************+')
disp(['List of selected points = ',num2str((setCandidates')),'  (for all clusters)'])
 WEIGHTS = zeros(length(setCandidates),length(ECMdata_cluster_test)) ; 
for icluster = 1:length(ECMdata_cluster_test)
    wDECM = ECMdata_cluster_test{icluster}.wRED; 
    zDEM = ECMdata_cluster_test{icluster}.setPoints;
    
    [II,JJ] = ismember(zDEM,setCandidates) ; 
    
end

%disp(num2str(setElements'))
%   clipboard('copy',num2str(setElements'));