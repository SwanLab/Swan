function [ECMdata_cluster,setCandidates,setElements] = ECMhyperLOCALperMODE(BstRED_l,BasisPone,DATA,wSTs,DATAoffline)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
% JAHO, 2-Jan-2024, Balmes 185, Barcelona 
% ---------------------------------
if nargin  == 0
    load('tmp.mat')
end


nmodes = size(BstRED_l,2) ;


disp('SEQUENTIAL METHOD FOR HYPERREDUCTION')
disp('********************************************************************')
disp(['Determining basis matrices for internal forces for each mode'])
disp('********************************************************************')
BasisFint_permode = cell(1,nmodes);

for imode = 1:nmodes
     
    
    BasisFint_permode{imode} = QbasisMatrixIntegrand(BstRED_l(:,imode),BasisPone,DATA,wSTs,DATAoffline) ;
    disp(['Disp.  mode = ',num2str(imode),';  N modes FINT =  ',num2str(size(BasisFint_permode{imode},2))]) ;
end

% Sort the number of columns of each BasisFint_permode
DATAoffline = DefaultField(DATAoffline,'OrderSequentialECM','SORTED') ; % If =2,
[dummy,nmodesFINT] = cellfun(@size,BasisFint_permode) ;

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
        
        TOTAL_number_perm =factorial(length(BasisFint_permode)); 
        
        if  TOTAL_number_perm <= DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM
        NEW_ORDER_clusters = perms(1:length(BasisFint_permode)) ; 
        else
        NEW_ORDER_clusters = zeros(DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM,length(BasisFint_permode)) ;
        for  irandom = 1:DATAoffline.NUMBER_RANDOM_RUNS_OrderSequentialECM
            NEW_ORDER_clusters(irandom,:) = randperm(length(BasisFint_permode)) ;
        end
        
        end
    case 'RANDN_fixed'
        rng(1)
        NEW_ORDER_clusters = randperm(length(BasisFint_permode)) ;
        %  NEW_ORDER_clusters = 1:length(BasisFint_permode) ;
    case 'PRINCIPAL_ANGLES_SUBSPACES'
        % Criterion based on the principla angles formed by the column
        % spaces of BasisFint_permode
       NEW_ORDER_clusters  =  SortECM_sequentialANGLES(BasisFint_permode) ; 
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
    
    [ECMdata_cluster_test,setCandidates_test] = ECM_local_clusteringTAILORED(nmodes,NEW_ORDER_clusters(itrials,:),wSTs,DATAoffline,BasisFint_permode,...
        DATA)  ;
    
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


disp(['Total number of points (had they been  selected independently) = ',num2str(sum(nmodesFINT))])

disp(['Total number of ECM candidate points = ',num2str(length(setCandidates))])


setElements = large2smallREP(setCandidates,DATA.MESH.ngaus) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements (for all modes)'])

disp(num2str(setElements'))
%   clipboard('copy',num2str(setElements'));