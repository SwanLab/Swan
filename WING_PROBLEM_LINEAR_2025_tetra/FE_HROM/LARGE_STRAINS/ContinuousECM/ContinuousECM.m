function ECMdata= ContinuousECM(BstRED_l,BasisPone,DATA,wSTs,DATAoffline,DATA_GENGAUSS,...
    MESH,Nst)
% JAHO, 13-JAN-2021
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/01_2D_beam_LARGE/README_OECM.pdf
if nargin == 0
    load('tmp.mat')
end
ECMdata = [] ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'DIRECT_POLYNOMIAL_FITTING_INTEGRAND',1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matrix of reduced internal forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNAPredFINT_nw = BasisF_from_BasisStress_PK1(BstRED_l,BasisPone,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;

% This is for storing the information for the CECM (to be analyzed outside)
DATAoffline =DefaultField(DATAoffline,'SAVE_CECM_INFO',0); % = 1;
    ndim = size(MESH.COOR,2) ;

NAME_LOC_WSDECM = 'DATA_DECM.mat' ;
if DATAoffline.SAVE_CECM_INFO == 1
    A = SNAPredFINT_nw' ;  
    W = wSTs;    
      
    COOR = MESH.COOR' ;
    xINI = Nst*COOR(:);
    xINI = reshape(xINI,ndim,[])';
    xLIM = zeros(ndim,2) ;
    for idim = 1:ndim
        xLIM(idim,1) = min(MESH.COOR(:,idim)) ;
        xLIM(idim,2) = max(MESH.COOR(:,idim)) ;
    end
    save(NAME_LOC_WSDECM,'A','W','xINI','xLIM')
    clear A
end

SNAPredFINT = bsxfun(@times,SNAPredFINT_nw,sqrt_wST) ;
% Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
%DATAoffline.errorFINT = 1e-3;
DATAsvd.RELATIVE_SVD = 1;
[Q,S_FINT,V_FINT,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;


Vw = bsxfun(@times,V_FINT',S_FINT)' ; 
IntApproxSVD = (Q*Vw')'*sqrt(wSTs) ; 
INTexac = SNAPredFINT_nw'*wSTs ; 
ErrorINtAPPROXsvd = norm(INTexac-IntApproxSVD)/norm(INTexac)*100 ; 
disp(['Error associated to the SVD =',num2str(ErrorINtAPPROXsvd),' % (for a SVD tolerance =',num2str(DATAoffline.errorFINT)  ,' )'])






% % Enlarge the basis matris for SNAPredFINT
a  = sqrt_wST - Q*(Q'*sqrt_wST) ;
TOL  = 1e-10 ;
if norm(a) > TOL
    INCLUDE_ADDITIONAL_COLUMN = 1;
else
    INCLUDE_ADDITIONAL_COLUMN = 0 ;
end

%VOL_rel =  norm(a)/norm(sqrt_wST) ;
%disp(['Relative norm of the orthogonal complement  of sqrt(W) to Q  =',num2str(VOL_rel),'  (if it is close to 1, it is mandatory to augment the matrix)']) ;
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'ENFORCE_SUM_WEIGHTS_EQUAL_VOLUME',1);
if  DATA_GENGAUSS.ENFORCE_SUM_WEIGHTS_EQUAL_VOLUME == 1 && INCLUDE_ADDITIONAL_COLUMN ==1
    a = a/norm(a) ;
    Q = [a,Q] ;   % It takes less in converge....
    %  Q = [sqrt_wST,Q] ;
    VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT = 3;
else
    VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT = 0 ;
end


Vw = Q'*SNAPredFINT ; 
IntApproxSVD = (Q*Vw)'*sqrt(wSTs) ; 
INTexac = SNAPredFINT_nw'*wSTs ; 
ErrorINtAPPROXsvd = norm(INTexac-IntApproxSVD)/norm(INTexac)*100 ; 
disp(['Error associated to the SVD, enriched with the constant function =',num2str(ErrorINtAPPROXsvd),' % (for a SVD tolerance =',num2str(DATAoffline.errorFINT)  ,' )'])

DATAoffline = DefaultField(DATAoffline,'USE_SELECTIVE_DECM',1) ; 
DATAoffline = DefaultField(DATAoffline,'ListElementsExclude_fromGID',[]) ;

DATA_GENGAUSS.USE_SELECTIVE_DECM =  DATAoffline.USE_SELECTIVE_DECM;
DATA_GENGAUSS.ListElementsExclude_fromGID =DATAoffline.ListElementsExclude_fromGID  ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POINTS CANDIDATE FOR THE DISCRETE ECM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[INDSEL,DATAOUTdistorsion] = PointsToIncludeECM(DATA_GENGAUSS,MESH,Nst,DATA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Empirical cubature method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_ECM = [] ;
DATA_ECM.TOL = DATAoffline.errorECM ;
DATA_ECM.IND_POINTS_CANDIDATES = INDSEL ;

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'USE_SELECTIVE_DECM',1) ; % = 1;
if DATA_GENGAUSS.USE_SELECTIVE_DECM == 0
    % Version before 3-DEc-2021
    [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;
else
    % New version (after 3-DEc-2021)
    [setPoints,wRED,~,~]= DiscreteEmpiricalCubatureMethod(Q',wSTs,DATA_ECM)  ;
    
end

IntegralExact =SNAPredFINT_nw'*wSTs  ;
IntegrationError = IntegralExact - SNAPredFINT_nw(setPoints,:)'*wRED ;
IntegrationError = norm(IntegrationError)/norm(IntegralExact)*100;
disp(['Actual integration error using DECM= ',num2str(IntegrationError),' % (prescribed tolerance fint =',num2str(DATAoffline.errorFINT*100), '%'])



% Determining the indices of the associated elements
ngausTOTAL = length(wSTs);
nstrain = size(BstRED_l,1)/ngausTOTAL;
ngaus = ngausTOTAL/size(MESH.CN,1) ;
setElements = large2smallREP(setPoints,ngaus) ;
disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
%clipboard('copy',num2str(setElements'));

HYPERREDUCED_VARIABLES.setPoints = setPoints ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = wRED ;  % Set associated WEights
HYPERREDUCED_VARIABLES.nstrain = nstrain ;
HYPERREDUCED_VARIABLES.PHI = bsxfun(@times,Q,1./sqrt(wSTs)) ;
HYPERREDUCED_VARIABLES.nstrain = nstrain ;

% Printing reduced set of elements
%--------------------------------------
DATA_GENGAUSS = ECMpointsPRINT(DATA_GENGAUSS,MESH,HYPERREDUCED_VARIABLES) ;
% -------------------------------------

% ------------------------------------------
% OTHER INPUTS FOR THE CONTINUOUS ALGORITHM
%-------------------------------------------
if DATA_GENGAUSS.DIRECT_POLYNOMIAL_FITTING_INTEGRAND ==0
VAR_SMOOTH_FE.BdomRED_gauss = BstRED_l;
VAR_SMOOTH_FE.BasisS_gauss = BasisPone;
end
VAR_SMOOTH_FE.ngausE = size(DATA.MESH.posgp,2);
% Inverse of the singular values of BasisF times its right-singular vectors
% --------------------------------------------------------------------------
if DATA_GENGAUSS.DIRECT_POLYNOMIAL_FITTING_INTEGRAND ==0
VAR_SMOOTH_FE.invSVsingular_F = bsxfun(@times,V_FINT',1./S_FINT)' ;
end
% Matrix of connectivities
% ---------------------------
VAR_SMOOTH_FE.COOR = MESH.COOR ;
% Matrix of coordinates
% ---------------------------
VAR_SMOOTH_FE.CN = MESH.CN ;
% Type of element
% --------------------------
VAR_SMOOTH_FE.TypeElement = MESH.TypeElement ;
% Dela. triangulation
% --------------------------
VAR_SMOOTH_FE.DELTRIANG = delaunayTriangulation(MESH.COOR);
% ORDER POLYNOMIALS
% ----------------------------
nnodeE = size(MESH.CN,2) ;
% 

VAR_SMOOTH_FE.DIRECT_POLYNOMIAL_FITTING_INTEGRAND = DATA_GENGAUSS.DIRECT_POLYNOMIAL_FITTING_INTEGRAND ; 
 

DATA = DefaultField(DATA,'posgp_given',[]) ; 
 
ngausSTRESS = size(DATA.MESH.posgp,2) ;
norder_poly = round(ngausSTRESS^(1/ndim)-1); 
ORDER_POLYNOMIALS = norder_poly*ones(1,ndim);
switch VAR_SMOOTH_FE.TypeElement
    case 'Quadrilateral'
        IND_POLYG_ELEMENT = [1 2 3 4 1] ;   
         %   ORDER_POLYNOMIALS  =[norder_poly,norder_poly] ;     
    case 'Hexahedra'
        IND_POLYG_ELEMENT = [1 2 3 4 5 6 7 8 1] ;
    case 'Triangle'
        IND_POLYG_ELEMENT = [1 2 3 1] ;         
    otherwise
        error('element not implemented')
end
VAR_SMOOTH_FE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
VAR_SMOOTH_FE.IND_POLYG_ELEMENT = IND_POLYG_ELEMENT;   % Local numbering of corner nodes (polygon)
VAR_SMOOTH_FE.nstrain = HYPERREDUCED_VARIABLES.nstrain  ;
VAR_SMOOTH_FE.wSTs =  wSTs ;

EXACTint = HYPERREDUCED_VARIABLES.PHI'*wSTs ;
APPROXint = HYPERREDUCED_VARIABLES.PHI(HYPERREDUCED_VARIABLES.setPoints,:)'*HYPERREDUCED_VARIABLES.WdomRED ;
eee = norm(EXACTint-APPROXint)/norm(EXACTint)*100 ;
disp(['DECM ***> Difference between approx. and ex. int  of LEFT SINGULAR VECTORS =',num2str(eee),' %']);


DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_INTERNAL_FORCE_MODES',0) ;

if DATA_GENGAUSS.PLOT_INTERNAL_FORCE_MODES == 1
    DATALOCfint.NameFileMesh = DATA_GENGAUSS.NameFileMesh_FINT ;
    DATALOCfint.MaterialType= MESH.MaterialType ;
    GidPostProcessModesFINT(MESH.COOR,MESH.CN,...
        MESH.TypeElement,HYPERREDUCED_VARIABLES.PHI,DATA.MESH.posgp,DATALOCfint);
end

VAR_SMOOTH_FE.BasisS_nodes = [] ;


%DATAoffline =DefaultField(DATAoffline,'SAVE_CECM_INFO',1); % = 1;
NAME_LOC_WSDECM = 'DATA_DECM.mat' ;
if DATAoffline.SAVE_CECM_INFO == 2
    A = SNAPredFINT_nw' ;  
    W = wSTs;   
    PHI = HYPERREDUCED_VARIABLES.PHI ; 
    wDECM = HYPERREDUCED_VARIABLES.WdomRED;
    zDECM =   HYPERREDUCED_VARIABLES.setPoints ;
    ndim = size(MESH.COOR,2) ;
    COOR = MESH.COOR' ;
    xINI = Nst*COOR(:);
    xINI = reshape(xINI,ndim,[])';
    xLIM = zeros(ndim,2) ;
    for idim = 1:ndim
        xLIM(idim,1) = min(MESH.COOR(:,idim)) ;
        xLIM(idim,2) = max(MESH.COOR(:,idim)) ;
    end
    save(NAME_LOC_WSDECM,'A','W','xINI','xLIM','wDECM','zDECM','PHI')
    clear A
end
 


DATAIN = [] ;
[xNEW,wNEW,DATAIN,ELEMENTS_xNEW,HISTORY_ITERATIONS,VAR_SMOOTH_FE,Ninterpolation] =  GeneralizedGaussLARGE...
    (wSTs,MESH,DATAIN,HYPERREDUCED_VARIABLES,Nst,DATA_GENGAUSS,VAR_SMOOTH_FE)  ;


%%%
ECMdata.COOR = xNEW ;
ECMdata.wRED = wNEW ;
ECMdata.setPoints = []  ;
ECMdata.setElements = ELEMENTS_xNEW ;
ECMdata.Ninterpolation = Ninterpolation ;







