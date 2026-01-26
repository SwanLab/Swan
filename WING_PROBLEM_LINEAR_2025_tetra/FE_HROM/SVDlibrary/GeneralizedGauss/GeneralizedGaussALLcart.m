function [xGAUSS,wGAUSS,DATAapprox] = GeneralizedGaussALLcart(G,W,xMAT,z,w,DATALOC,S,V)
% Generalized gaussian cubature
%dbstop('6')
if nargin ==0
    load('tmp2.mat')
end

DATALOC = DefaultField(DATALOC,'PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS',0) ;
DATALOC = DefaultField(DATALOC,'TOLERANCE_FORCE_POINTS_TO_REMAIN_WITHIN_THE_DOMAIN',0) ;
DATALOC = DefaultField(DATALOC,'SAVE_XGAUSS',1 ) ;

xINT = [] ;
% MATRIX of original basis functions
% -------------------------------------------------------------------------
G = bsxfun(@times,G',1./sqrt(W))' ;
PHI_der = []  ;
% Exact integral
% ---------------
b= PHI'*W ;
% Starting Cubature rule
% ----------------------
wOLD =w ; xOLD = xMAT(z,:) ; DATALOC.mINI = length(z) ;
DATALOC.zOLD =  z; POINTS_all.x = cell( length(z),1) ; 
POINTS_all.w = cell( length(z),1) ; POINTS_all.z = cell( length(z),1) ;

%VSinv = bsxfun(@times,V',1./S)' ;
ALLPOSITIVE = [] ;
CONVERGENCE = 1 ;
iter = 1;
while CONVERGENCE ==1 && length(xOLD)>1    
    [xNEW,wNEW,CONVERGENCE,DATALOC] = Iter_RemovePointALLDcart(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC,...
        VSinv) ;
    if CONVERGENCE == 1
        xOLD = xNEW ;
        wOLD = wNEW ;
        DATALOC.zOLD = DATALOC.zNEW ;
        POINTS_all.x{iter} = xNEW ;
        POINTS_all.w{iter} = wNEW ;
        POINTS_all.z{iter} =  DATALOC.zNEW ;
        ALLPOSITIVE(end+1) = all(wNEW>0) ;
        iter = iter + 1 ;
    end
end

if ~isempty(ALLPOSITIVE)
    IND_ALL = find(ALLPOSITIVE==1) ;
    IND_LAST = IND_ALL(end) ;
    xOLD = POINTS_all.x{IND_LAST} ;
    wOLD = POINTS_all.w{IND_LAST} ;
    zOLD = POINTS_all.z{IND_LAST} ;
    POINTS_all.x =  POINTS_all.x(1:IND_LAST) ;
    POINTS_all.w =  POINTS_all.w(1:IND_LAST) ;
    POINTS_all.z =  POINTS_all.z(1:IND_LAST) ;
end


disp('--------------------------------------------------')
disp(['Final integration rule with m =',num2str(length(xOLD)),' POINTS  (of ',num2str(length(z)),'). ',...
    '. Rank Basis = ',num2str(length(b))])
disp('--------------------------------------------------')

disp('Integration error')

PHIk_y = EvaluateBasisFunctionAnalytical(xOLD,VSinv,DATALOC,0,1) ;


bNEW = PHIk_y'*wOLD ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error (%) =',num2str(errorINT)])

%---------------
xINITIAL =xMAT(z,:) ;
wINITIAL = w ;
wGAUSS = wOLD ;
xGAUSS = xOLD ;
if DATALOC.SAVE_XGAUSS ==1
    save(DATALOC.NAMEWS,'xINITIAL','wINITIAL','wGAUSS','xGAUSS','errorINT','POINTS_all','DATALOC');
end


ndim = size(xOLD,2) ;

if  ndim ==3
    Plot3DpointsCECM(DATALOC,xINITIAL,xOLD,xMAT,wOLD,wINITIAL,zOLD) ;
else
    error('Option not implemnted')
end

DATAapprox = [] ;