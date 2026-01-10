function [PHI,COORsort,Wnew,DATAIN,xx,INDinterpolate] = ...
    OrderGaussPointsNonCartesian(DATALOC,COORg,z,PHI_OLD,W,w,DATAIN,COOR,CN,TypeElement,xxOLD)
%  COUNTERPART OF OrderGaussPointsCartesian FOR NON CARTESIAN MESHES
% -----------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
Wnew=[];
% -------------------------------------------
% STEP 1: CONSTRUCTING THE CARTESIAN MESH
% ------------------------------------------
% How to decice the size of the cartesian mesh ???? It is dictated by
% an input variable, and by the coordinates of the ECM points

xx = CartesianMeshForInterpolation(COORg,DATALOC,z,xxOLD) ; 


[xINT{1} xINT{2}]= meshgrid(xx{1},xx{2})  ;  % Matrices xINT and yINT for interpolation
COORsort = [xINT{1}(:),xINT{2}(:)];  % Coordinates ---> Sorted in a cartesian fashion

%---------------------------------------
% STEP 2 --> Indices ECM points
%------------------------------------------
% NO LONGER NEEDED
% TriangECM = delaunayTriangulation(COORsort(:,1),COORsort(:,2));  %
% [zNEW,distanceNEW ]= nearestNeighbor(TriangECM, COORg(z,1),COORg(z,2));  % Nearest points old coordinates
% if sum(distanceNEW ~=0)
%     error('Distances must be zero !!!! The ECM points does  not coincide with any point of the discretization')
% end

 
% STEP 2:  Initialization
% --------------------------------------------------------------------------
PHI =zeros(size(COORsort,1),size(PHI_OLD,2)) ;    % New matrix of modes
% --------------------------------------------------------------------------
% STEP 3: Determine which points of COORsort are outside the domain
%---------------------------------------------------------------------------
% FOR THIS, WE NEED TO DEFINE THE BOUNDARY OF OUR DOMAIN AS A POLYGON.  
%  
COORpolygonBND = PolygonBoundary2DfromGID(DATALOC) ; 
% Next we determine which points are outside the domain. 
figure(456)
hold on
plot(COORpolygonBND(:,1),COORpolygonBND(:,2),'LineWidth',2) ; 
axis equal
plot(COORsort(:,1),COORsort(:,2),'r.')
[INDinside,INDon] = inpolygon(COORsort(:,1),COORsort(:,2),COORpolygonBND(:,1),COORpolygonBND(:,2)) ; 

% INTERPOLATION SHALL BE CARRIED OUT ONLY ON  the indexes INDinside
% ------------------------------------------------------------------
INDinterpolate = 1:size(COORsort,1) ; 
INDinterpolate  = INDinterpolate(INDinside) ; 

% STEP 4:  INTERPOLATION 
% ------------------------
% Type of interpolation 
%-------------------------
%DATALOC.InterpolationTypeNonCartesianMesh = 'Quadratic'; 'Linear' ;  
DATALOC = DefaultField(DATALOC,'InterpolationTypeNonCartesianMesh','Quadratic') ;  

switch DATALOC.InterpolationTypeNonCartesianMesh 
    case 'Quadratic'
        InterpolationMethod = 'natural' ; 
    case 'Linear'
        InterpolationMethod = 'linear' ; 
        case 'Nearest'
              InterpolationMethod = 'nearest' ; 
    otherwise
        error('Option not implemented')
end
ExtrapMethod = 'linear' ; % In principle, it is not needed

for imode = 1:size(PHI,2)     
    if imode == 1 
        Finterp = scatteredInterpolant(COORg(:,1),COORg(:,2),PHI_OLD(:,imode),InterpolationMethod,ExtrapMethod) ;  
    else
        Finterp.Values = PHI_OLD(:,imode) ; 
    end
    
    PHI(INDinterpolate,imode) = Finterp(COORsort(INDinterpolate,1),COORsort(INDinterpolate,2)) ; 
end

disp('---')
 
 
%
%     %------------------------------------
%     % Therefore
%     % PHI --> Matrix with the functions we wish to integrate
%     PHI =  PHI_OLD(INDnear,:) ; %
%     % Check that everything is consistent
%     % Exact integral
%     EXACT = PHI'*Wnew;
%     % Approx. integral
%     APPROX = PHI(zNEW,:)'*w;
%     ERROR_chech = norm(EXACT-APPROX)/norm(EXACT)*100
%     DATAIN.MSGPRINT{end+1}  = ['Int. Error after changing indexes = ',num2str(ERROR_chech),'%'] ;
%     disp(DATAIN.MSGPRINT{end})