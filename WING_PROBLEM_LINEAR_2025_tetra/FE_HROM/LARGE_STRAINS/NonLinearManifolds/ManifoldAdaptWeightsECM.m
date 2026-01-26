function [wALL,setPointsALL,qLATENT,SNAPfint] = ManifoldAdaptWeightsECM(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT)
% JAHO, Manifold-adaptive weights ECM, 16th September 2025, UPC, TErrassa
% % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
if nargin == 0
    load('tmp.mat')
    DATAoffline.IntegrateExactlyLinearResponse_SAW_ECM = 0 ;
    DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM = [] ;
    close all
    format long g
end
%
%  methodology AS DESCRIBED IN
%   See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
% The functions we wish to integrate can be represented as a sequence of vectors of R^M, M being the number of Gauss points:
% f_1, f_2....f_n, where f_i = f(q_i), q_i < q_i+1
% The vectors are normalized with respect to the vector of FE weights:  f_i^T*W*f_i = 1.
% From a geometric point of view, the above sequence can be interpreted as (smooth ) rotations in R^N of the initial vector f_1.
% Suppose that up to certain  1 < i < n, the ECM returns a common point z_1 for f_1, f_2 ...f_i.
% Then, when it moves to q_{i+1}, the algorithm is  unable to find a positive weight such that f_{i+1}(z_1)*w_{i+1} = f_{i+1}^T W.
% It is at this juncture that the "new" ingredient comes into play:
% rather than simply enlarge the set of candidate points in order to find a integration point with positive weight,
% we intend to keep the previous point z_1, and simply seek another point z_2 such that w_1 and w_2 positive.
%Yet this can be only be done if we increase the size of the basis functions of cluster i+1.
% The simplest idea in this respect is to form the new basis with f_{i+1} and the orthogonal complement of f_i (with
% respect to f_{i+1}).
%

% We remove the case qLATENT = 0  (on the grounds that the modes here may be pure noise)
[qLATENT,ii] = sort(qLATENT) ; % We sort qLATENT in ascending order
[aab,bbb] = find(abs(qLATENT)<1e-10) ;
if ~isempty(bbb)
    qLATENT(bbb) = [] ;
    ii(bbb) = [] ;
end
SNAPfint = SNAPfint(:,ii) ;  % Integrand snapshot cell array (each cell = one variable of qLATENT)

sqW = sqrt(OPERFE.wSTs) ; % Square root FE Gauss weights
sqWfe = diag(sparse(sqW)) ; % Diagonal, sparse matrix containing the FE Gauss weights (square root)
DATAoffline = DefaultField(DATAoffline,'IntegrateExactlyLinearResponse_SAW_ECM',0) ;

ModesFint_default_w = {} ;
TolerancesFintDEFAULT_w = [] ;

if  DATAoffline.IntegrateExactlyLinearResponse_SAW_ECM == 1
    % INTEGRATING EXACTLY THE LINEAR RESPONSE
    % HOW TO DETERMINE A "NOISELESS" LINEAR RESPONSE
    NSAMPLES_LINEAR = min(3,length(SNAPfint)); % hyperparameter, number of samples
    % in which the system is considered linear
    [~,indexLINEAR] = sort(abs(qLATENT));
    indexLINEAR=  indexLINEAR(1:NSAMPLES_LINEAR) ;
    SNAPfint_LINEAR=  SNAPfint(indexLINEAR) ;
    nCOMP = size(SNAPfint_LINEAR{1},2) ;
    MODES_LINEAR_FINT = zeros(size(SNAPfint_LINEAR{1})) ;
    SNAPfint_LINEAR = cell2mat(SNAPfint_LINEAR) ;
    for icomp  = 1:nCOMP
        compSELECT= icomp:nCOMP:size(SNAPfint_LINEAR,2) ;
        [UUU,SSS,VVV] =  SVDT(SNAPfint_LINEAR(:,compSELECT));
        MODES_LINEAR_FINT(:,icomp) = UUU(:,1) ;
    end
    
    MODES_LINEAR_FINT = SVDT(sqWfe*MODES_LINEAR_FINT)  ;
    
    ModesFint_default_w{end+1} =  MODES_LINEAR_FINT ;
    TolerancesFintDEFAULT_w(end+1)  = 0 ;
    
else
    MODES_LINEAR_FINT = [] ;
end


DATAoffline = DefaultField(DATAoffline,'IntegrateExactlyVolume_SAW_ECM',1) ;
DATAoffline = DefaultField(DATAoffline,'Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM',[]) ;

SNAPfint_w  = cell(1,size(SNAPfint,2)) ;

NumberModesCluster = zeros(size(SNAPfint_w)) ;
INCLUDE_sum_vol  =DATAoffline.IntegrateExactlyVolume_SAW_ECM;

if ~isempty(DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM)
    disp(['SVD entire internal force matrix (can be speed up by using direct randomization)'])
    [UUUdom,SSSdom,~] =   SRSVD(sqWfe*cell2mat(SNAPfint));
    nmodesINC = min(DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM,length(SSSdom)) ;
    ModesSVD_weighted = UUUdom(:,1:nmodesINC) ;
    ModesFint_default_w{end+1} =  ModesSVD_weighted ;
    TolerancesFintDEFAULT_w(end+1)  = 0 ;
end


%if DATAoffline.IntegrateExactlyVolume_SAW_ECM == 1
ModesFint_default_w{end+1} =  sqW ;
TolerancesFintDEFAULT_w(end+1)  = 0 ;
%end

% Ubasis_w_const = SRSVD(ModesFint_default_w,TolerancesFintDEFAULT_w) ;  % This is the "initial" basis function
SetCandidates=  1:length(sqW) ;
Ubasis_w_const = cell2mat(ModesFint_default_w) ;


DATA_ECM.TOL = 0 ;
setPoints_cluster = cell(length(SNAPfint),1) ;
wRED_cluster = cell(length(SNAPfint),1) ;
DATA_ECM.TOL =0;

DIR_qNEG = find(qLATENT<= 0) ;
DIR_qPOS = find(qLATENT> 0) ;


[aaa,ind_sort_neg] = sort(qLATENT(DIR_qNEG),'descend') ;
DIR_qNEG_sorted = DIR_qNEG(ind_sort_neg) ;

[aaa,ind_sort_pos] = sort(qLATENT(DIR_qPOS)) ;
DIR_qPOS_sorted = DIR_qPOS(ind_sort_pos) ;

ORDER_clusters = [DIR_qNEG_sorted,DIR_qPOS_sorted] ;


for iclusterLOC = 1:length(SNAPfint)
    
    icluster = ORDER_clusters(iclusterLOC) ;
    
    disp(['iclusterLOC = ',num2str(iclusterLOC),' GLO = ',num2str(icluster),' q = ',num2str(qLATENT(icluster))]) ; 
    
    
    
    SNAPlocI =  sqWfe*SNAPfint{icluster} ; % This is the sbapshot matrix of the current cluster
    % Next we apply the partitioned SVD, forcing to include the "previous"
    % basis matrix (by using a tolerance 0)
    TOL = [0,DATAoffline.errorFINT] ;
    DDD.HIDE_OUTPUT = 1;
    
    wRED = [] ;
    iterFAIL = 0 ;
    while  isempty(wRED)
        if iterFAIL >0
            TOLloc = zeros(size(TOL)) ;
        else
            TOLloc = TOL ;
        end
        [UU,SS,VV] = SRSVD({Ubasis_w_const,SNAPlocI},TOLloc,DDD) ;
        DATA_ECM.IND_POINTS_CANDIDATES =  SetCandidates ;
        [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_NoIter(UU',OPERFE.wSTs,DATA_ECM)  ;
        
        if isempty(wRED)
            % We have to enlarge the basis matrix
            Ubasis_w_const = Ubasis_w_prev ;
            iterFAIL = iterFAIL + 1;
        else
            
            wRED_cluster{icluster} = wRED ;
            setPoints_cluster{icluster}= setPoints ;
            %   if  icluster ==1
            SetCandidates =  setPoints ;
            %  else
            %     SetCandidates = unique([SetCandidates(:);setPoints(:)]);
                
            % end
            
            Ubasis_w_prev = UU ;
            
        end
        disp(['Total number of candidate points =',num2str(length(SetCandidates))])
    end
    
    NumberModesCluster(icluster) = length(SS) ;
end

figure(42)
hold on
xlabel('Snapshot')
ylabel('Number of internal force modes')
bar(NumberModesCluster)



%%%%%%%%%%%%%%%%%5
%Creating the chart qLATENT-wECM,zECM
setPointsALL = unique(cell2mat(setPoints_cluster)) ;

% for icluster = 1:length(ECMdata_cluster)
%     setPointsALL{icluster} =  ECMdata_cluster{icluster}.setPoints  ;
% end
% setPointsALL  =unique(cell2mat(setPointsALL' )) ;

ncluster = length(wRED_cluster) ;
npointsALL = length(setPointsALL) ;
wALL = zeros(npointsALL,ncluster) ;



for icluster = 1:ncluster
    setPloc = setPoints_cluster{icluster}  ;
    [dummy1,III,JJJ] = intersect(setPloc,setPointsALL,'stable') ;
    wALL(JJJ,icluster) = wRED_cluster{icluster}  ;
    
    
end

figure(4532)
hold on
title('WEIGHTS VERSUS LATENT VARIABLE (HYPERELASTIC PROBLEM)')
xlabel('q')
ylabel('w')
for iddd  = 1:size(wALL,1)
    eeee = large2smallINCLUDEREP(setPointsALL(iddd),DATA.MESH.ngaus) ;
    plot(qLATENT,wALL(iddd,:),'DisplayName',['Point =',num2str(setPointsALL(iddd)),' Elem = ',num2str(eeee)])
end
legend show
