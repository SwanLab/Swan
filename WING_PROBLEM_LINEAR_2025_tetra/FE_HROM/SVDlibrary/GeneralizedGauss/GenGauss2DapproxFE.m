function [xGAUSS,wGAUSS,DATALOC] = GenGauss2DapproxFE(PHI,W,COORg,xx,z,w,DATALOC,coorECM,VAR_SMOOTH_FE)
% Generalized gaussian cubature
% See Generalized Gaussian quadrature rules on arbitrary polygons
%dbstop('6')
if nargin ==0
    load('tmp1.mat')
    DATALOC.RemoveInaddmissiblePointsAtTheEnd = 0 ;
end

% DATALOC = DefaultField(DATALOC,'LOAD_PREVIOUSLY_COMPUTED_POINTS',0) ;
% DATALOC = DefaultField(DATALOC,'CRITERION_ELIMINATION',0) ;
% DATALOC = DefaultField(DATALOC,'TOL_low',DATALOC.TOL)  ;
% DATALOC = DefaultField(DATALOC,'NAMEWS_diary','DIARY.txt')  ;
% DATALOC = DefaultField(DATALOC,'xZ',[])  ;
% DATALOC = DefaultField(DATALOC,'FUN_GRADIENT_DIRECTLY',[])  ;
%
% DATALOC = DefaultField(DATALOC,'SHOW_FIGURES',1 ) ;
% DATALOC = DefaultField(DATALOC,'SAVE_XGAUSS',1 ) ;


%diary([DATALOC.NAMEWS_diary])

%DATALOC = DefaultField(DATALOC,'NEGATIVE_check_during_iterations',1) ;

% Spatial gridW
% ------------------------------
if isempty(VAR_SMOOTH_FE)
    %  Methods based on an underlying cartersian GRID
    xINT = cell(1,2) ;
    [xINT{1}, xINT{2}]= meshgrid(xx{1},xx{2})  ;  % Matrices xINT and yINT for interpolation
else
    xINT=[] ;
end

% MATRIX of original basis functions (assuming that J = Lambda'*sqrt(W))
% -------------------------------------------------------------------------
% %dbstop('62')
% %if   DATALOC.PHIisTRANSPOSE ==0
%     PHI = bsxfun(@times,J',1./sqrt(W)) ;
%     clear J
% %else
% %    PHI = J' ;
% %end


DATALOC = DefaultField(DATALOC,'APPROX_FUN__DERI',[]) ;
%DATALOC.APPROX_FUN__DERI = DefaultField(DATALOC.APPROX_FUN__DERI,'ACTIVE',1) ; % Replace by 2D fitting
DATALOC.APPROX_FUN__DERI.PHI = PHI ;

DATALOC.APPROX_FUN__DERI= DefaultField(DATALOC.APPROX_FUN__DERI,'METHOD','LOCAL_FITTING') ;  % or SVD_based_FITTING
DATALOC.APPROX_FUN__DERI= DefaultField(DATALOC.APPROX_FUN__DERI,'SVD_INTERPOLATION',[]) ;



switch DATALOC.APPROX_FUN__DERI.METHOD
    case 'SVD_based_FITTING'
        % METHOD BASED ON THE SVD DECOMPOSITION OF EACH COLUMN OF PHI
        % Method based on interpolation via SVD
        [DATAfitSVD ] = SVDinterpolation_FITTING(DATALOC.APPROX_FUN__DERI.SVD_INTERPOLATION,PHI,xx{1},xx{2}) ;
        DATALOC.APPROX_FUN__DERI.SVD_INTERPOLATION.DATAfitSVD = DATAfitSVD;
        
        DATALOC.APPROX_FUN__DERI.COOR = [] ;
        % case 'FE_INTERPOLATION'
        
        
    otherwise
        DATALOC.APPROX_FUN__DERI.COOR = COORg ;
        
end
DATALOC =  DefaultField(DATALOC,'MSGPRINT',{}) ;
%DATAapprox = DATALOC.APPROX_FUN__DERI ;
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
if isempty(VAR_SMOOTH_FE)
    hx = xx{1}(2)-xx{1}(1);
    hy = xx{2}(2)-xx{2}(1) ;
    DATALOC.hx = hx ; DATALOC.hy = hy ;
    DATALOC.xLIM(1,:) =[xx{1}(1),xx{1}(end)] ;
    DATALOC.xLIM(2,:) =[xx{2}(1),xx{2}(end)] ;
end



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
% - --------------

b=   DATALOC.ExactIntegral; % PHI'*W ;


% Starting Cubature rule
% ----------------------
wOLD =w ;
%if isempty(DATALOC.xZ)
%if isempty(z)
xOLD = coorECM ;
 
DATALOC.PATHiniPOINTS = cell(size(xOLD,1),size(xOLD,1)) ;

 

%else
%    xOLD = COORg(z,:) ;
%end

% else
%     xOLD = DATALOC.xZ ;
% end
DATALOC.mINI = size(xOLD,1) ;

POINTS_all.x = {} ;
POINTS_all.w = {} ;
DATALOC.dxRELmax =0 ;
DATALOC.dyRELmax =0 ;

%VSinv = bsxfun(@times,V',1./S)' ;

% DATALOC = DefaultField(DATALOC,'RemoveInaddmissiblePointsAtTheEnd',0) ;
%
% if DATALOC.RemoveInaddmissiblePointsAtTheEnd == 1
%     % Created 8-Apr2020, to test the possibility of removing the inadd.
%     % points at the end of the iteration procedure
%     [xOLD,wOLD,DATALOC,POINTS_all,iter] = ...
%         Iter_RemPnt_AtTheEnd(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC) ;
%
%
%
% else
CONVERGENCE = 1 ;
iter = 1;
DATALOC.iter = iter;

for ipoint = 1:size(xOLD,1)
    DATALOC.PATHiniPOINTS{ipoint,iter } = xOLD(ipoint,:) ;
end

INDEX_xNEW_respect_XOLD =  1:size(xOLD,1) ; 
% At each iteration of the ensuing loop, a point of xOLD is removed 
% We shall asign to the set of initial points the labels 1,2....npoints 
% At each iteration, we have to ascertain which are the points that remains
% in the set 
if ~isempty(VAR_SMOOTH_FE)
[~,~,POLYINFOloc]=     EvaluateBasisFunctionAtX_FEinterp(xOLD,[],VAR_SMOOTH_FE,[]) ; 
ELEMENTS_TO_PLOT =POLYINFOloc.ELEMENTS_CONTAINING_xNEW  ; 

end

while CONVERGENCE ==1 & length(xOLD)>1
    
    
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO] = Iter_RemPnt2Dapprox(xOLD,wOLD,xINT,PHI,PHI_der,b,DATALOC,VAR_SMOOTH_FE) ;
    if CONVERGENCE == 1
       
       if ~isempty(POLYINFO)
           ELEMENTS_TO_PLOT = [ELEMENTS_TO_PLOT;POLYINFO.ELEMENTS_CONTAINING_xNEW  ] ; 
       end
     %   POINTS_all.x{iter} = xNEW ;
     %   POINTS_all.w{iter} = wNEW ;
     iter = iter + 1 ;
     DATALOC.iter = iter ;
     INDEX_xNEW_respect_XOLD(DATALOC.REMOVED_INDEX) = [] ;
     for  ipoint = 1:size(xNEW,1)
         DATALOC.PATHiniPOINTS{INDEX_xNEW_respect_XOLD(ipoint),iter } = xNEW(ipoint,:) ;
     end
     xOLD = xNEW ;
     wOLD = wNEW ;
    end
end

%end




DATALOC.MSGPRINT{end+1}='--------------------------------------------------' ;
disp(DATALOC.MSGPRINT{end})
DATALOC.MSGPRINT{end+1} = ['Final integration rule with m =',num2str(length(xOLD)),' POINTS  (of ',num2str(size(coorECM,1)),'). ',...
    '. Rank Basis = ',num2str(length(b))];
disp(DATALOC.MSGPRINT{end})
DATALOC.MSGPRINT{end+1}='--------------------------------------------------' ;
disp(DATALOC.MSGPRINT{end})


DATALOC.MSGPRINT{end+1}='Integration error' ;

disp(DATALOC.MSGPRINT{end})
%     if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%         PHIk_y = zeros(length(wOLD),size(PHI,2));
%         for i=1:size(PHI,2)
%             PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%             PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xOLD(:,1),xOLD(:,2),'cubic') ;
%         end
%     else
%PHIk_y = EvaluateBasisFunctionAtX(xOLD,VSinv,DATALOC,0,1) ;
[PHIk_y,~, POLYINFO]= EvaluateBasisFunctionAtX_approx(xOLD, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

%   end

bNEW = PHIk_y'*wOLD ;
errorINT = norm(bNEW-b)/norm(b)*100;
DATALOC.MSGPRINT{end+1}=['Error (%) =',num2str(errorINT)] ;
disp(DATALOC.MSGPRINT{end}) ;

%---------------
xINITIAL =coorECM; %(z,:) ;
wINITIAL = w ;
wGAUSS = wOLD ;
xGAUSS = xOLD ;
% if DATALOC.SAVE_XGAUSS ==1
%     save(DATALOC.NAMEWS,'xINITIAL','wINITIAL','wGAUSS','xGAUSS','errorINT','POINTS_all','DATALOC')
% end

DATALOC.SHOW_FIGURES = 1;
if DATALOC.SHOW_FIGURES==1
    
    
    
    figure(15)
    hold on
    xlabel('x')
    ylabel('y')
    %  if DATALOC.INITIAL_DISCRETIZATION_TENSORPRODUCT ==1
    
    if  ~isempty(POLYINFO)
        % Plot Elements affected during the iterations
        % ---------------------------------------------
%         INDused = 1:length(POLYINFO.COEFFSpolynomial) ;
%         INDEMPT = cellfun(@isempty,POLYINFO.COEFFSpolynomial) ;
%         INDused(INDEMPT) = [];
        INDused= unique(ELEMENTS_TO_PLOT) ; 
        for ielem = 1:length(INDused)
            elemLOC = INDused(ielem) ;
            NODESloc  = VAR_SMOOTH_FE.CN(elemLOC,VAR_SMOOTH_FE.IND_POLYG_ELEMENT)' ;
            COORloc = VAR_SMOOTH_FE.COOR(NODESloc,:) ;
            plot(COORloc(:,1),COORloc(:,2),'LineWidth',1) ;
        end
        
        % Plotting path during removal iterations 
        % DATALOC.PATHiniPOINTS
        
        for ipoint =1:size(DATALOC.PATHiniPOINTS,1) 
            locPATH = DATALOC.PATHiniPOINTS(ipoint,:) ; 
            indN = cellfun(@isempty,locPATH) ; 
            locPATH = locPATH(indN==0) ; 
            locPATH = cellfun(@transpose,locPATH,'UniformOutput',false) ;  
            COOR = cell2mat(locPATH) ; 
 %           vxlabels = arrayfun(@(n) {sprintf('k%d', n)}, (1:size(COOR,2))');
%Hpl = text(COOR(1,:), COOR(2,:), vxlabels,'HorizontalAlignment','center');
 
            
            plot(COOR(1,:),COOR(2,:),'Marker','.') ; 
            
            
            
        end
        
        
        
        
    else
        DATALOC = DefaultField(DATALOC,'INDinterpolate',[]) ;
        if ~isempty(DATALOC.INDinterpolate)
            plot(COORg(DATALOC.INDinterpolate,1),COORg(DATALOC.INDinterpolate,2),'.','Color',[0 1 0]) ;
        else
            plot(COORg(:,1),COORg(:,2),'.','Color',[0 1 0]) ;
        end
    end
    
    %  end
    h = plot(xINITIAL(:,1),xINITIAL(:,2),'o','Color',[0 0 1],'MarkerSize',8);
    h2 = plot(xOLD(:,1),xOLD(:,2),'x','Color',[1 0 0],'MarkerSize',8);
    LLL{1} = ['Before opt., m =',num2str(size(xINITIAL,1))] ;
    LLL{2} = ['After opt., m =',num2str(size(xOLD,1))] ;
    legend([h,h2],LLL);
    axis equal
    
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
dt = DelaunayTri(COORg(:,1),COORg(:,2));
NP = nearestNeighbor(dt, xOLD(:,1),xOLD(:,2));

xNEAR = COORg(NP,:); % Nearest point
% % Integration error
% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     PHIk_y = zeros(size(xNEAR,1),size(PHI,2));
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%         PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEAR(:,1),xNEAR(:,2),'cubic') ;
%     end
% else
%PHIk_y = EvaluateBasisFunctionAtX(xNEAR,VSinv,DATALOC,0,1) ;

PHIk_y = EvaluateBasisFunctionAtX_approx(xNEAR, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

%end


% Weights
wNEAR = PHIk_y'\b ;
bNEAR = PHIk_y'*wNEAR ;
errorINTnear = norm(bNEAR-b)/norm(b)*100;

DATALOC.MSGPRINT{end+1}=['Error nearest points (%) =',num2str(errorINTnear)] ;
disp(DATALOC.MSGPRINT{end}) ;
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
