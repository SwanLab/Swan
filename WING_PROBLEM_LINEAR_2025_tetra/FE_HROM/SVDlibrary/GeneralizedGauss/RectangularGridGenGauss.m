function [xx,coorECM,DATALOC,DATAIN,PHI,zNEW,COORsort,Wnew] =...
    RectangularGridGenGauss(w,COORg,HYPERREDUCED_VARIABLES,W,z,DATALOC,DATAIN,COOR,CN,TypeElement)

xx = cell(1,2) ;
TOLloc = 1e-6; % TOLERANCE for determining which points are aligned with those  with y =  ymin

idim = 1 ; iOTHER = 2;
xmin = min(COORg(:,idim)) ;
DIFFx = abs(COORg(:,idim)-xmin) ;
IROW = find(DIFFx <= TOLloc);
xx{iOTHER} = sort(COORg(IROW,iOTHER)) ;

idim = 2 ; iOTHER = 1;
xmin = min(COORg(:,idim)) ;
DIFFx = abs(COORg(:,idim)-xmin) ;
IROW = find(DIFFx <= TOLloc);
xx{iOTHER} = sort(COORg(IROW,iOTHER)) ;
DATALOC.ExactIntegral = HYPERREDUCED_VARIABLES.PHI'*W;
INDinterpolate = [] ;

coorECM = COORg(z,:) ; % COordinates initial set of points

%  DATALOC = DefaultField(DATALOC,'CENTRAL_GAUSS_POINTS_ONLY_IN_Q9',0) ;
% if DATALOC.CENTRAL_GAUSS_POINTS_ONLY_IN_Q9 ==1
%     % Only consider central gauss points in Q9 elements
%    switch TypeElement
%        case 'Quadrilateral'
%            if  size(CN,2) == 9
%                % Q9 ELEMENTS
%                nelem = size(CN,1) ;
%                ngausELEM = 9 ;
%                GaussPointsAll = small2large(1:nelem,ngausELEM) ;
%                GaussSelect = GaussPointsAll(ngausELEM:ngausELEM:end) ;
%                COORg = COORg(GaussSelect,:) ;
%                HYPERREDUCED_VARIABLES.PHI = HYPERREDUCED_VARIABLES.PHI(GaussSelect,:) ;
%            end
%    end
%
% end

DATALOC=DefaultField(DATALOC,'USE_CARTESIAN_GENERAL_METHOD',0);
DATALOC=DefaultField(DATALOC,'APPROX_FUN__DERI',[]);
DATALOC.APPROX_FUN__DERI=DefaultField(DATALOC.APPROX_FUN__DERI,'METHOD','SVD_based_FITTING');


switch DATALOC.APPROX_FUN__DERI.METHOD
    
    case {'FE_INTERPOLATION'}
        % No cartesian mesh is needed. The integrand and its derivatives
        % are computed via FE evaluation
        xx = [];   zNEW = z ; 
        PHI = HYPERREDUCED_VARIABLES.PHI ; 
        COORsort = COORg ; 
        Wnew = W ; 
         
    
    case {'SVD_based_FITTING','LOCAL_FITTING'}
        
        % We have to construct a cartesian grid for interpolation purposes 
    
    if length(xx{1})*length(xx{2}) ~= size(COORg,1)  || DATALOC.USE_CARTESIAN_GENERAL_METHOD ==1
        
        % The FE mesh is not cartesian ---or it is cartesian, yet we ignore
        % it by setting DATALOC.USE_CARTESIAN_GENERAL_METHOD =1
        DATAIN.MSGPRINT{end+1} ='NON-CARTESIAN MESH' ;
        disp(DATAIN.MSGPRINT{end})
        
        zNEW = [] ;
        
        [PHI,COORsort,Wnew,DATAIN,xx,INDinterpolate] = ...
            OrderGaussPointsNonCartesian(DATALOC,COORg,z,HYPERREDUCED_VARIABLES.PHI,W,w,...
            DATAIN,COOR,CN,TypeElement,xx) ;
    else
        DATAIN.MSGPRINT{end+1} ='CARTESIAN MESH' ;
        disp(DATAIN.MSGPRINT{end})
        % CARTESIAN MESH.
        % It is assumed
        %that the Gauss points are located in a cartesian mesh.
        %Therefore, all we have to do is to change the indexes of the coordinates
        %COORg, the ECM points, the matrix of snapshots PHI, and  the FE weights.
        
        [PHI,zNEW,COORsort,Wnew,DATAIN] = OrderGaussPointsCartesian(xx,COORg,z,HYPERREDUCED_VARIABLES.PHI,W,w,DATAIN) ;
        
    end
    
    DATALOC.INDinterpolate = INDinterpolate;
    
end

