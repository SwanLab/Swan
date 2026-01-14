function  [xINInew,xx,INDEXperm,MPOINTS] = ArrangePointsCartesianFromFEmesh2D(xINI)
if nargin == 0
    load('tmp1.mat')
end
ndim = 2;
xx = cell(1,2) ;
TOLloc = 1e-4; % TOLERANCE for determining which points are aligned with those  with y =  ymin
%  figure(1)
%  hold on
%  plot(xINI(:,1),xINI(:,2),'x')

%
idim = 1 ; iOTHER = 2;
xmin = min(xINI(:,idim)) ;
xmax  = max(xINI(:,idim)) ;
dx = xmax-xmin ;
DIFFx = abs(xINI(:,idim)-xmin)/dx ;
IROW = find(DIFFx <= TOLloc);
xx{iOTHER} = sort(xINI(IROW,iOTHER)) ;

idim = 2 ; iOTHER = 1;
xmin = min(xINI(:,idim)) ;
xmax  = max(xINI(:,idim)) ;
dx = xmax-xmin ;
DIFFx = abs(xINI(:,idim)-xmin)/dx ;
IROW = find(DIFFx <= TOLloc);
xx{iOTHER} = sort(xINI(IROW,iOTHER)) ;
[MPOINTS,~] = cellfun(@size,xx) ; 


[xINInew{1} xINInew{2}]= meshgrid(xx{1},xx{2})  ;  % Matrices xINT and yINT for interpolation
xINInew = [xINInew{1}(:),xINInew{2}(:)];  % Coordinates ---> Sorted in a cartesian fashion

 
% Now we seek a set of indices INDnear such that COORg(INDnear,:) =
% COORsort
dt = delaunayTriangulation(xINI(:,1),xINI(:,2));  % Del. triangulation old coordinates
INDEXperm = nearestNeighbor(dt, xINInew(:,1),xINInew(:,2));  % Nearest points new coordinates

EEE = norm(xINInew-xINI(INDEXperm,:),'fro') ; 
NNN = norm(xINI,'fro') ; 

ERROR_rel = EEE/NNN ; 
TOL = 1e-10 ; 
if ERROR_rel >TOL 
    error('Non-cartesian mesh')
end
 





%
%         % We have to construct a cartesian grid for interpolation purposes
%
%     if length(xx{1})*length(xx{2}) ~= size(xINI,1)  || DATALOC.USE_CARTESIAN_GENERAL_METHOD ==1
%
%         % The FE mesh is not cartesian ---or it is cartesian, yet we ignore
%         % it by setting DATALOC.USE_CARTESIAN_GENERAL_METHOD =1
%         DATAIN.MSGPRINT{end+1} ='NON-CARTESIAN MESH' ;
%         disp(DATAIN.MSGPRINT{end})
%
%         zNEW = [] ;
%
%         [PHI,COORsort,Wnew,DATAIN,xx,INDinterpolate] = ...
%             OrderGaussPointsNonCartesian(DATALOC,xINI,z,HYPERREDUCED_VARIABLES.PHI,W,w,...
%             DATAIN,COOR,CN,TypeElement,xx) ;
%     else
%         DATAIN.MSGPRINT{end+1} ='CARTESIAN MESH' ;
%         disp(DATAIN.MSGPRINT{end})
%         % CARTESIAN MESH.
%         % It is assumed
%         %that the Gauss points are located in a cartesian mesh.
%         %Therefore, all we have to do is to change the indexes of the coordinates
%         %xINI, the ECM points, the matrix of snapshots PHI, and  the FE weights.
%
%         [PHI,zNEW,COORsort,Wnew,DATAIN] = OrderGaussPointsCartesian(xx,xINI,z,HYPERREDUCED_VARIABLES.PHI,W,w,DATAIN) ;
%
%     end
%
%     DATALOC.INDinterpolate = INDinterpolate;
%
% end
%
