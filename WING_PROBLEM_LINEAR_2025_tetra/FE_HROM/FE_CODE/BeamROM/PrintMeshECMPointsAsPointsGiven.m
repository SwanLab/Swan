function [COOR,CNredASPOINTS]=  PrintMeshECMPointsAsPoints(COOR,CNred,DATAIN,DATA_REFMESH,setElements);

if nargin == 0
    load('tmp1.mat')
end

% S1) Finding the centroid of each ECM element. In principle, this requires
% having at one's disposal the shape functions of each node.
ndim = size(COOR,2) ;
%nelem = 1;
nnodeE = size(CNred,2) ;
%[ Nelem,posgp  ] = ComputeNelemALL(DATA_REFMESH.TypeElement,nnodeE,ndim,nelem) ; % Shape functin of a single element
%ngaus = size(posgp,2);
%W = diag(DATA_REFMESH.WdiagRHS) ;
%W = W(1:ndim:end) ;  % Weights

% COOR_cube =COOR(CNred(1,:),:) ; 
% centro = sum(COOR_cube,1)/size(COOR_cube,1) ; 
% COOR_cube = bsxfun(@minus,COOR_cube',centro')' ; 

switch  DATA_REFMESH.TypeElement
    case 'Hexahedra'
        if nnodeE == 27
            [~,~,~,~,COOR_cube] = Hexahedra27NInPoints ;
        else
            
            COOR_cube =[-1  1  1    -1   -1  1    1  -1
                -1  -1  1   1    -1  -1   1   1
                -1  -1  -1  -1   1   1    1   1]' ;
            
        end
    otherwise
        error('Option not implemented')
end
     
% EStimate lenght 
% ------------------
dL = norm(COOR(CNred(1,1),:)-COOR(CNred(1,2),:)) ;  
%COOR_cube = bsxfun(@minus,COOR_cube',centro')' ; 
DATAIN = DefaultField(DATAIN,'SCALE_FACTOR_FOR_PRINTING_ECM_POINTS',0.1) ; %

TOL = DATAIN.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS ; 
COOR_cube = COOR_cube*TOL*dL ; 
COORgiven = DATAIN.COORDINATES_POINTS_TO_PRINT ; 
nnodesGIVEN = size(COORgiven,1) ; 
nnodesADD = nnodesGIVEN*size(COOR_cube,1) ; 
 
%DATAIN = DefaultField(DATAIN,'LengthCube_ECMpointsPRINT',dL) ; %
%COOR_cube = COOR_cube*DATAIN.LengthCube_ECMpointsPRINT ;
ADDitionalCOOR = zeros(nnodesADD,ndim) ;
CNredASPOINTS = zeros(size(CNred)) ;

numNODES = size(COOR,1) ;
numNODESloc = 0 ;

for iecm = 1:size(CNred,1)
  %  CNloc = CNred(iecm,:) ;
    % DOFS = small2large(CNloc,ndim) ;
    % Computing centroid
    % ------------------
  %  iecmGLO =  setElements(iecm) ;  % Global element
    % Computing CENTROID
    % -------------------
  %  GAUSS = small2large(iecmGLO, ngaus);
  %  Wloc = W(GAUSS) ;
  %  VOL = sum(Wloc) ;
  %  COORnodes = COOR(CNloc,:)';
  %  CENTROIDall = Nelem*COORnodes(:);
  %  CENTROID = zeros(ndim,1) ;
  %  for idim = 1:ndim
  %      CENTROID(idim) = Wloc'*CENTROIDall(idim:ndim:end)/VOL ;
  %  end
    % -------------------------------------------------------------
    % Now we have to create an hexahedra that encloses the computed
    % CENTROID. Input length paramater --> DATAIN.LengthCube_ECMpointsPRINT
    % Creating new points
    for innode = 1:nnodeE
        numNODES = numNODES +1 ; % Label of the incoming node
        numNODESloc = numNODESloc +1 ;
        ADDitionalCOOR(numNODESloc,:) = COORgiven(iecm,:)+COOR_cube(innode,:) ;
        CNredASPOINTS(iecm,innode) = numNODES ;
    end
    
    
end

% if ~isempty(DATAIN.COORDINATES_POINTS_TO_PRINT)
%     ADDitionalCOOR = DATAIN.COORDINATES_POINTS_TO_PRINT ; 
% end
COOR = [COOR;ADDitionalCOOR ] ;
