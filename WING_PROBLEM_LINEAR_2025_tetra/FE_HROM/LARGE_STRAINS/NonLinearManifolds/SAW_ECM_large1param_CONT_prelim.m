function [ECMdata_cluster,setCandidates,qLATENT,SNAPfint] = SAW_ECM_large1param_CONT_prelim(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT)
% JAHO, 13-Sept-2025, Saturday, Thriller, Caffe bar, Sarajevo, Bosnia.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
if nargin == 0
    load('tmp.mat')
    DATAoffline.IntegrateExactlyLinearResponse_SAW_ECM = 0; 
    DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM = 1; 
    close all
end


%SNAPfint = cell2mat(SNAPfint) ;
% We remove the case qLATENT = 0  (on the grounds that the modes here may be pure noise)
[qLATENT,ii] = sort(qLATENT) ;
[aab,bbb] = find(abs(qLATENT)<1e-10) ;
if ~isempty(bbb)
    qLATENT(bbb) = [] ;
    ii(bbb) = [] ;
end
SNAPfint = SNAPfint(:,ii) ;

sqW = sqrt(OPERFE.wSTs) ;
sqWfe = diag(sparse(sqW)) ;
DATAoffline = DefaultField(DATAoffline,'IntegrateExactlyLinearResponse_SAW_ECM',0) ;

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
end


for icluster = 1:size(SNAPfint,2)
    
    disp(['icluster = ',num2str(icluster)])
    %     if icluster == 299
    %         disp('')
    %     end
    %     icluster_back = max(1,icluster-DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    %     icluster_forw = min(size(SNAPfint,2),icluster+DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    %     icluster_select = [icluster_back:icluster_forw] ;
    
    SNAPlocI =  SNAPfint{icluster} ;
    % INCLUDE ELASTIC SNAPSHOTS/integ. exact volum
    SNAPlocI = sqWfe*SNAPlocI ;
    if isempty(MODES_LINEAR_FINT) && isempty(DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM)
        % Linear modes are not included
        if INCLUDE_sum_vol ==0
            if  size(SNAPlocI,2) >1
                TOL = [DATAoffline.errorFINT] ;
                DDD.ISRELATIVE = 1;
                [UU,SS,VV] = SVDT(SNAPlocI,TOL,DDD) ;
            else
                UU = SNAPlocI/norm(SNAPlocI,'fro') ;
                SS = 1;
            end
        else
            TOL = [0,DATAoffline.errorFINT] ;
            DDD.HIDE_OUTPUT = 1;
            [UU,SS,VV] = SRSVD({sqW,SNAPlocI},TOL,DDD) ;
        end
    elseif  ~isempty(DATAoffline.Number_FirstSVDmodesFINT_includeLOCALBASES_SAW_ECM)
         TOL = [0,0,DATAoffline.errorFINT] ;
        DDD.HIDE_OUTPUT = 1;
        
        [UU,SS,VV] = SRSVD({sqW,ModesSVD_weighted,SNAPlocI},TOL,DDD) ;
    else
        % Linear modes are included . In this option, we also enforce
        % constant mode
        %     if INCLUDE_sum_vol ==0
        %          TOL = [0,DATAoffline.errorFINT] ;
        %  DDD.HIDE_OUTPUT = 1;
        
        % [UU,SS,VV] = SRSVD({MODES_LINEAR_FINT,SNAPlocI},TOL,DDD) ;
        %  else
        TOL = [0,0,DATAoffline.errorFINT] ;
        DDD.HIDE_OUTPUT = 1;
        
        [UU,SS,VV] = SRSVD({sqW,MODES_LINEAR_FINT,SNAPlocI},TOL,DDD) ;
        % end
    end
    
    SNAPfint_w{icluster} =  UU ;
    NumberModesCluster(icluster) = length(SS) ;
end

figure(42)
hold on
xlabel('Snapshot')
ylabel('Number of internal force modes')
bar(NumberModesCluster)


%
[ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ...
    SAW_ECM(OPERFE.wSTs,DATAoffline,SNAPfint_w,DATA) ;


%
% warning('BORRAR ESTO, ESTOY ENGAÑANDO AL PROGRAMA')
%     load('tmpWORKecm.mat','setPoints','wRED')
% ECMdata.setPoints = setPoints ;
% ECMdata.wRED.Values = repmat(wRED,1,length(qLATENT)) ;