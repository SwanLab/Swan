function [ECMdata,DATAOUT] = ECMcluster_global(BasisFint_cluster,...
    OPERFE,DATA,DATAoffline) ; 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/PAPER_SAW_ECM.mlx
if nargin == 0
    load('tmp.mat')
end
%
wSTs = OPERFE.wSTs ; 
sqrt_wST = sqrt(wSTs) ;
for i = 1:length(BasisFint_cluster)
    BasisFint_cluster{i} = bsxfun(@times, BasisFint_cluster{i},sqrt_wST) ; 
end
 
  DATAlocc = [] ; 

DATAoffline = DefaultField(DATAoffline,'TOLERANCE_GLOBAL_FINT_FOR_GLOBAL_CLUSTERING',1e-8)  ; 
TOL = DATAoffline.TOLERANCE_GLOBAL_FINT_FOR_GLOBAL_CLUSTERING ; 

[Q,S,V] = SRSVD(BasisFint_cluster,TOL,DATAlocc) ; 

disp(['Number of basis vectors of Q = ',num2str(size(Q,2))]) ;    

% % Enlarge the basis matris for SNAPredFINT
a  = sqrt_wST - Q*(Q'*sqrt_wST) ;
if norm(a) > 1e-10
    a = a/norm(a) ;
    Q = [Q,a] ;
%    SNAPredFINT_nw = [ones(size(a)),SNAPredFINT_nw] ; %
    
end


TIMEdecm = tic ; 

DATAoffline = DefaultField(DATAoffline,'ListElementsExclude_fromGID',[]) ; 
if ~isempty(DATAoffline.ListElementsExclude_fromGID)
     ListElementsToExclude = load(DATAoffline.ListElementsExclude_fromGID) ; 
    ListElementsToExclude = ListElementsToExclude(:,1) ; 
    ngausELEM = DATA.MESH.ngaus_STRESS; 
    ListGaussToExclude = small2large(ListElementsToExclude,ngausELEM) ; 
    INDSEL  = setdiff(1:size(Q,1),ListGaussToExclude) ; 
else
    INDSEL  = 1:size(Q,1) ; 
end

DATA_ECM.IND_POINTS_CANDIDATES = INDSEL ;
DATA_ECM.TOL = DATAoffline.errorECM ;
DATAoffline = DefaultField(DATAoffline,'USE_SELECTIVE_DECM',1) ; % = 1;
if DATAoffline.USE_SELECTIVE_DECM == 0
    % Version before 3-DEc-2021
    [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;
else
    % New version (after 3-DEc-2021)
    [setPoints,wRED,ERROR_GLO,DATAOUT]= DiscreteEmpiricalCubatureMethod(Q',wSTs,DATA_ECM)  ;
    
end

DATAOUT.TIMEdecm = toc(TIMEdecm) ; 

disp(['Time to select the integration points =',num2str(DATAOUT.TIMEdecm)])


    ECMdata.setPoints = setPoints ;
    ECMdata.wRED = wRED ;
%     proporP = length(setPoints)/DATA.MESH.ngausT*100;
%     disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
    setElements = large2smallREP(setPoints,DATA.MESH.ngaus) ;
%     disp('****************************+')
     disp(['List of selected m_e = ',num2str(length(setElements)),' elements'])
     disp(num2str(setElements'))
%    % clipboard('copy',num2str(setElements'));
    ECMdata.setElements = setElements ;
 
    
    