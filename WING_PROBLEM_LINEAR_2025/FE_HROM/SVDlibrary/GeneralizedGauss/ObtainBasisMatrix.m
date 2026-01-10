function [Lambda,S,V,INTexac,DATA] = ObtainBasisMatrix(Xf,W,DATA)
% SVD of Xf*sqrt(W)
% --------------------------------------
%dbstop('5')
if nargin == 0
    load('tmp1.mat')
end
% Exact integral
% ---------------
if ~iscell(Xf)
    INTexac = Xf'*W ; % High-fidelity integral of all snapshots
else
    INTexac = cell(size(Xf)) ;
    for imat = 1:length(Xf)
        if  ischar(Xf{1})
            load(Xf{imat},'Aloc') ;
            INTexac{imat} =Aloc'*W ;
        else
            INTexac{imat} =Xf{imat}'*W ;
        end
    end    
    INTexac = cell2mat(INTexac') ;
end

V = sum(W) ;    % Volume of the domain

if ~iscell(Xf)
    Xf = bsxfun(@times,Xf,sqrt(W)) ;
else    
    for ipart = 1:length(Xf)
        Xf{ipart} = bsxfun(@times,Xf{ipart},sqrt(W)) ;
    end    
end


% Orthogonal basis matrix Lambda for the column space of Xf 
disp('SVD of Aw ')
 
DATA = DefaultField(DATA,'UsePartitionedRandomizedAlgorithm',1) ;

DATA = DefaultField(DATA,'PointsToBeIntegratedWithHighAccuracy',[]) ; 

if ~isempty(DATA.PointsToBeIntegratedWithHighAccuracy)
    [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    SegretateSVDdecm(Xf, DATA ) ; 
else
    if  ~iscell(Xf)
        if DATA.UsePartitionedRandomizedAlgorithm == 1
            DATAsvd.EPSILON_GLO = DATA.TOLsvdXf ;
            [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    RSVDqp(Xf, DATAsvd.EPSILON_GLO ) ;
        else
            DATAsvd.RELATIVE_SVD =1 ;
            [Lambda,S,V] =    SVDT(Xf,DATA.TOLsvdXf,DATAsvd ) ;
        end
    else
        if isempty(DATA.epsilon)
            epsilon = 0*ones(size(Xf)) ;
        else
            epsilon = DATA.epsilon ;
        end
        [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    RSVDqp(Xf,epsilon,DATAsvd) ;
    end
end


% % Enlarge the basis matris for SNAPredFINT
a  = sqrt(W) - Lambda*(Lambda'*sqrt(W)) ;
TOL  = 1e-16 ;
if norm(a) > TOL
    INCLUDE_ADDITIONAL_COLUMN = 1;
else
    INCLUDE_ADDITIONAL_COLUMN = 0 ;
end
DATA = DefaultField(DATA,'Evaluate_Fun_Gradients_via_FITTING',0) ;
if    INCLUDE_ADDITIONAL_COLUMN ==1 && DATA.Evaluate_Fun_Gradients_via_FITTING==1
    % We have not implemented for the case
    % DATA.Evaluate_Fun_Gradients_via_FITTING==0 (analytical evaluation).
    % In this case, it should be enforced by including the constant
    % function in the integrand
    a = a/norm(a) ;
    Lambda = [a,Lambda] ;   % It takes less in converge....
    %  Q = [sqrt_wST,Q] ;
    %     S = S/S(1)  ;
    %     S = [1; S] ;
    
end




Vw = Lambda'*Xf ; 
IntApproxSVD = (Lambda*Vw)'*sqrt(W) ; 
ErrorINtAPPROXsvd = norm(INTexac-IntApproxSVD)/norm(INTexac)*100 ; 
disp(['Error associated to the SVD =',num2str(ErrorINtAPPROXsvd),' % (for a SVD tolerance =',num2str(DATA.TOLsvdXf)  ,' )'])

% SVD residual 
% ----------------------
 
EXAMINE_RESIDUALS = 0;
if EXAMINE_RESIDUALS == 1
    E_svd = Xf - Lambda*Vw ;
    % Residual for each point
    % --------------------------
    E_svd = sum(E_svd.^2,2) ;
    % Relative residual
    Xf = sum(Xf.^2,2) ;
    relRESID = sqrt(E_svd./Xf) ;
    DATA.ExcludePointsPercentageHighSVDResidual = 1 ;
    nPOINTS_exclude =  20 ;  ceil(length(relRESID)*DATA.ExcludePointsPercentageHighSVDResidual/100) ;
    [aa,bb] = sort(relRESID) ;
    ExcludedPoints = bb(1:nPOINTS_exclude) ;
    
    % load('tmp3.mat','INDEXperm')
    %     ExcludedPoints_old  = INDEXperm(ExcludedPoints) ;
    %    eOLD = large2small(ExcludedPoints_old,9)
    %     clipboard('copy',num2str(eOLD'))
    
    DATA.ECM_POINTS_INCLUDE = setdiff(DATA.ECM_POINTS_INCLUDE,ExcludedPoints) ;
    
    
end



PLOT_SVDerror = 0; 
if  PLOT_SVDerror ==1 
figure(356)
hold on
title('SVD error ')
SingVsq =  (S.*S) ;
SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
%normEf2 = sqrt(sum(SingVsq) - cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
normEf2 = sqrt(cumsum(SingVsq)) ;
normEf2 = sort(normEf2,'descend');
hold on
lS = log10(normEf2);

h = plot([0:length(S)-1],lS,'r');

xlabel('Number of Modes')
ylabel('log(error)')

%legend(h)
end

    
     

DATA = DefaultField(DATA,'TruncationCriterionSVDmodesJumpOrderMagnitudeLOG',5) ;
if ~isempty(DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG)
    Sincre = log10(S(1:end-1)./S(2:end)) ;    
    III = find(Sincre >= DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG) ;    
    if ~isempty(III)
        nmodes= III(1) ;
        Lambda = Lambda(:,1:nmodes) ;
        S = S(1:nmodes) ;
        V = V(:,1:nmodes) ;
    end   
end

%%


 
DATAOUT.SingV_Jmin = S ;
DATAOUT.V_Jmin = V ;
DATAOUT.INTexac = INTexac;
