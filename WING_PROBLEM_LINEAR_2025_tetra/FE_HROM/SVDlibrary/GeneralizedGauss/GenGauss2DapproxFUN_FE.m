function [xGAUSS,wGAUSS,DATAapprox] = GenGauss2DapproxFUN_FE(J,W,xMAT,z,w,DATALOC,MPOINTS,S,V,XY)
% Generalized gaussian cubature
% See Generalized Gaussian quadrature rules on arbitrary polygons
%dbstop('6')
if nargin ==0
    
    load('tmp.mat')
    porder = DATALOC.PORDER ;
    
    %  load('tmp1.mat')
    
    
    
    
end

DATALOC = DefaultField(DATALOC,'LOAD_PREVIOUSLY_COMPUTED_POINTS',0) ;
DATALOC = DefaultField(DATALOC,'CRITERION_ELIMINATION',0) ;
DATALOC = DefaultField(DATALOC,'TOL_low',DATALOC.TOL)  ;
DATALOC = DefaultField(DATALOC,'NAMEWS_diary','DIARY.txt')  ;
DATALOC = DefaultField(DATALOC,'xZ',[])  ;
DATALOC = DefaultField(DATALOC,'FUN_GRADIENT_DIRECTLY',[])  ;

DATALOC = DefaultField(DATALOC,'SHOW_FIGURES',1 ) ;
DATALOC = DefaultField(DATALOC,'SAVE_XGAUSS',1 ) ;


%diary([DATALOC.NAMEWS_diary])

%DATALOC = DefaultField(DATALOC,'NEGATIVE_check_during_iterations',1) ;

% Spatial gridW
% ------------------------------
xx =   reshape(xMAT(:,1),MPOINTS(2),[]) ;
xx = xx(1,:)' ;
yy =  xMAT(1:MPOINTS(2),2) ;
[xINT{1} xINT{2}]= meshgrid(xx,yy)  ;  % Matrices xINT and yINT for interpolation

% MATRIX of original basis functions (assuming that J = Lambda'*sqrt(W))
% -------------------------------------------------------------------------
%dbstop('62')
%if   DATALOC.PHIisTRANSPOSE ==0
    PHI = bsxfun(@times,J',1./sqrt(W)) ;
    clear J
%else
%    PHI = J' ;
%end


DATALOC = DefaultField(DATALOC,'APPROX_FUN__DERI',[]) ;
DATALOC.APPROX_FUN__DERI = DefaultField(DATALOC.APPROX_FUN__DERI,'ACTIVE',0) ; % Replace by 2D fitting
DATALOC.APPROX_FUN__DERI.PHI = PHI ;

DATALOC.APPROX_FUN__DERI= DefaultField(DATALOC.APPROX_FUN__DERI,'METHOD','LOCAL_FITTING') ;  % or SVD_based_FITTING
DATALOC.APPROX_FUN__DERI= DefaultField(DATALOC.APPROX_FUN__DERI,'SVD_INTERPOLATION',[]) ;



switch DATALOC.APPROX_FUN__DERI.METHOD
    case 'SVD_based_FITTING'
% METHOD BASED ON THE SVD DECOMPOSITION OF EACH COLUMN OF PHI 
% Method based on interpolation via SVD
        [DATAfitSVD ] = SVDinterpolation_FITTING(DATALOC.APPROX_FUN__DERI.SVD_INTERPOLATION,PHI,xx,yy) ;
   DATALOC.APPROX_FUN__DERI.SVD_INTERPOLATION.DATAfitSVD = DATAfitSVD; 
 
      
        
end

DATAapprox = DATALOC.APPROX_FUN__DERI ; 
 % --------------------------------------------------------------


% Gradient of basis functions (numerical gradient)
% -------------------------------------------------------------------
% if  METHOD_der ==1
%     for idim = 1:length(XfDER)
%         PHI_der{idim} = XfDER{idim}*V ;
%         PHI_der{idim} = bsxfun(@times,PHI_der{idim}',(1./S))' ;
%         PHI_der{idim} = bsxfun(@times,PHI_der{idim},sqrt(W)) ;
%     end
% else
 hx = xx(2)-xx(1);
 hy = yy(2)-yy(1) ;
 DATALOC.hx = hx ; DATALOC.hy = hy ;
% % Computing numerical gradient
% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     
%     PHI_der{1} = zeros(size(PHI)) ;
%     PHI_der{2} = zeros(size(PHI)) ;
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),MPOINTS(2),[]);
%         [dx  dy]  =gradient(PHIloc,hx,hy) ;
%         
%         PHI_der{1}(:,i) = dx(:) ;
%         PHI_der{2}(:,i) = dy(:) ;
%         
%     end
%     
% else
    PHI_der = []  ;
%end
%
% end

%% Exact integral
% ---------------
b= PHI'*W ;

% Starting Cubature rule
% ----------------------
wOLD =w ;
if isempty(DATALOC.xZ)
    xOLD = xMAT(z,:) ;
else
    xOLD = DATALOC.xZ ;
end
DATALOC.mINI = length(z) ;

POINTS_all.x = {} ;
POINTS_all.w = {} ;
DATALOC.dxRELmax =0 ;
DATALOC.dyRELmax =0 ;

%VSinv = bsxfun(@times,V',1./S)' ;
if   DATALOC.LOAD_PREVIOUSLY_COMPUTED_POINTS ==0
    CONVERGENCE = 1 ;
    iter = 1;
    
    
    while CONVERGENCE ==1 & length(xOLD)>1
        [xNEW,wNEW,CONVERGENCE,DATALOC] = Iter_RemPnt2Dapprox(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC) ;
        if CONVERGENCE == 1
            xOLD = xNEW ;
            wOLD = wNEW ;
            POINTS_all.x{iter} = xNEW ;
            POINTS_all.w{iter} = wNEW ;
            iter = iter + 1 ;
        end
    end
    disp('--------------------------------------------------')
    disp(['Final integration rule with m =',num2str(length(xOLD)),' POINTS  (of ',num2str(length(z)),'). ',...
        '. Rank Basis = ',num2str(length(b))])
    disp('--------------------------------------------------')
    
    disp('Integration error')
%     if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%         PHIk_y = zeros(length(wOLD),size(PHI,2));
%         for i=1:size(PHI,2)
%             PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%             PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xOLD(:,1),xOLD(:,2),'cubic') ;
%         end
%     else
        %PHIk_y = EvaluateBasisFunctionAtX(xOLD,VSinv,DATALOC,0,1) ;
         PHIk_y = EvaluateBasisFunctionAtX_approx(xOLD, DATALOC.APPROX_FUN__DERI)  ;
        
 %   end
    
    bNEW = PHIk_y'*wOLD ;
    errorINT = norm(bNEW-b)/norm(b)*100;
    disp(['Error (%) =',num2str(errorINT)])
    
    %---------------
    xINITIAL =xMAT(z,:) ;
    wINITIAL = w ;
    wGAUSS = wOLD ;
    xGAUSS = xOLD ;
    if DATALOC.SAVE_XGAUSS ==1
        save(DATALOC.NAMEWS,'xINITIAL','wINITIAL','wGAUSS','xGAUSS','errorINT','POINTS_all','DATALOC')
    end
else
    load(DATALOC.NAMEWS,'xINITIAL','wINITIAL','wGAUSS','xGAUSS','errorINT','POINTS_all','DATALOC')
    disp(['Error (%) =',num2str(errorINT)])
    
end

if DATALOC.SHOW_FIGURES==1
    figure(15)
    hold on
    xlabel('x')
    ylabel('y')
    if DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT ==1
        plot(xMAT(:,1),xMAT(:,2),'.','Color',[0 1 0])
    end
    h = plot(xINITIAL(:,1),xINITIAL(:,2),'o','Color',[0 0 1],'MarkerSize',8);
    h2 = plot(xOLD(:,1),xOLD(:,2),'x','Color',[1 0 0],'MarkerSize',8);
    LLL{1} = ['Before opt., m =',num2str(size(xINITIAL,1))] ;
    LLL{2} = ['After opt., m =',num2str(size(xOLD,1))] ;
    legend([h,h2],LLL);
    
end

%dbstop('162')
% DATALOC = DefaultField(DATALOC,'HOLE',[]) ;
% 
% if ~isempty(DATALOC.HOLE)
%     switch  DATALOC.HOLE.TYPE
%         case 'POLYGONAL'
%             %             POINTS_pol = DATA.DATAREMOVEPOINTS.HOLE.POINTS ;
%             %             POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
%             %             INPol = inpolygon(x,y,POINTS_pol(:,1),POINTS_pol(:,2)) ;
%             %             Xf(INPol,:) = 0 ;
%             POINTS_pol = DATALOC.HOLE.POINTS ;
%             POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
%             for ipoints = 1:size(POINTS_pol,1)-1
%                 plot(POINTS_pol(ipoints:ipoints+1,1),POINTS_pol(ipoints:ipoints+1,2),'b')
%             end
%             
%         otherwise
%             error('Option not implemented')
%     end
% end





% Nearest points
% ----------------
dt = DelaunayTri(xMAT(:,1),xMAT(:,2));
NP = nearestNeighbor(dt, xOLD(:,1),xOLD(:,2));

xNEAR = xMAT(NP,:); % Nearest point
% % Integration error
% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     PHIk_y = zeros(size(xNEAR,1),size(PHI,2));
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%         PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEAR(:,1),xNEAR(:,2),'cubic') ;
%     end
% else
    %PHIk_y = EvaluateBasisFunctionAtX(xNEAR,VSinv,DATALOC,0,1) ;
    
     PHIk_y = EvaluateBasisFunctionAtX_approx(xNEAR, DATALOC.APPROX_FUN__DERI)  ;
    
%end


% Weights
wNEAR = PHIk_y'\b ;
bNEAR = PHIk_y'*wNEAR ;
errorINTnear = norm(bNEAR-b)/norm(b)*100;
disp(['Error nearest points (%) =',num2str(errorINTnear)])
%%%

if DATALOC.SHOW_FIGURES==1
    figure(16)
    hold on
    bar(sort(wOLD))
    ylabel('Weights (after optimization)')
    figure(17)
    hold on
    bar(sort(wINITIAL))
    ylabel('Weights (before optimization)')
    
end
%diary off
