function [xCECM,wCECM,zCECM,DATAapprox] = RemoveRecursivePointsCART(G,W,xINI,z,w,DATALOC,S,V)
% Remove points recursively from the initial set until arrive at the
% integration rule {xCECM,wCECM}
% Cartesian domains
if nargin ==0
    load('tmp1.mat')
end

DATALOC = DefaultField(DATALOC,'PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS',0) ;
DATALOC = DefaultField(DATALOC,'TOLERANCE_FORCE_POINTS_TO_REMAIN_WITHIN_THE_DOMAIN',0) ; % = 1;
DATALOC = DefaultField(DATALOC,'SAVE_XGAUSS',1 ) ;
% MATRIX of original basis functions
% -------------------------------------------------------------------------
G = bsxfun(@times,G',1./sqrt(W))' ;
% Exact integral
% ---------------
b= G*W ;

DATALOC = DefaultField(DATALOC,'Evaluate_Fun_Gradients_via_FITTING',0) ;
if DATALOC.Evaluate_Fun_Gradients_via_FITTING == 1
    DATAFITTING = GetCoefficientesFitting(xINI,G,DATALOC) ;
else
    DATAFITTING = [] ;
end
DATAapprox.DATAFITTING = DATAFITTING ;

% Starting Cubature rule
% ----------------------
wOLD =w ; xOLD = xINI(z,:) ; DATALOC.mINI = length(z) ;
DATALOC.zOLD =  z; POINTS_all.x = cell( length(z),1) ;
POINTS_all.w = cell( length(z),1) ; POINTS_all.z = cell( length(z),1) ;

VSinv = bsxfun(@times,V',1./S)' ;
ALLPOSITIVE = zeros(length(z),1) ;
CONVERGENCE = 1 ;
iter = 1;
while CONVERGENCE ==1 && size(xOLD,1)>1
    [xNEW,wNEW,CONVERGENCE,DATALOC] = Remove1PointCART(xOLD,wOLD,b,DATALOC,...
        VSinv,DATAFITTING) ;
    if CONVERGENCE == 1
        xOLD = xNEW ;
        wOLD = wNEW ;
        DATALOC.zOLD = DATALOC.zNEW ;
        POINTS_all.x{iter} = xNEW ;
        POINTS_all.w{iter} = wNEW ;
        POINTS_all.z{iter} =  DATALOC.zNEW ;
        ALLPOSITIVE(iter) = all(wNEW>0) ;
        iter = iter + 1 ;
    end
end

IND_ALL = find(ALLPOSITIVE==1) ;
if isempty(IND_ALL)
    xCECM = xOLD ;
    wCECM = wOLD ;
    zCECM = z;
else
    IND_LAST = IND_ALL(end) ;
    xCECM = POINTS_all.x{IND_LAST} ;
    wCECM = POINTS_all.w{IND_LAST} ;
    zCECM = POINTS_all.z{IND_LAST} ;
    POINTS_all.x =  POINTS_all.x(1:IND_LAST) ;
    POINTS_all.w =  POINTS_all.w(1:IND_LAST) ;
    POINTS_all.z =  POINTS_all.z(1:IND_LAST) ;
end



disp('--------------------------------------------------')
disp(['Final integration rule with m =',num2str(length(wCECM)),' POINTS  (of ',num2str(length(z)),')'])
disp('--------------------------------------------------')
disp('Integration error')
disp('-------------------------------------')

EVALUATE_GRADIENT = 0 ;
[gCECM,~ ]= Evaluate_Basis_Grad_Analytic(xCECM,VSinv,DATALOC,EVALUATE_GRADIENT,DATAFITTING)  ;


bNEW = gCECM*wCECM ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error (%) =',num2str(errorINT)])

%---------------
if DATALOC.SAVE_XGAUSS ==1
    save(DATALOC.NAMEWS,'wCECM','xCECM','errorINT','POINTS_all','DATALOC');
end

DATAapprox.INFO_iterations = POINTS_all;
DATAapprox.VSinv =VSinv ;
