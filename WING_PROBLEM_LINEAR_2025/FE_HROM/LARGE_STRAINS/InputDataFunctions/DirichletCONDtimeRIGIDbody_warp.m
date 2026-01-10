function      [U,a,DOFrLOC] =    DirichletCONDtimeRIGIDbody_warp(LOCVAR,DATA,ndim,MESH,GEOproperties) ;
% DIRICHLET BOUNDARY CONDITIONS,  beam type
% Seee /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
%
% iloadstate = 1;  % Loading stage
% LOCVAR.PRESCRIBED_DISP = [] ;
 
% LOCVAR.PRESCRIBED_DISP(iloadstate).AMPLITUDE = [disp_x,disp_y,disp_z,rot_x,rot_y,rot_z,warp] ;  %  
% LOCVAR.PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution

% 2D 
% LOCVAR.PRESCRIBED_DISP(iloadstate).AMPLITUDE = [disp_x,disp_y,rot_z] ;  %  
% LOCVAR.PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
% 23-March-2024, Barcelona, Balmes 185
%----------------------------------------------------
%
if nargin == 0
    load('tmp1.mat')
end

%%%%
isurface= LOCVAR.NUMBER_SURFACE ;  % Index whose surface is to be rotated + warped
rnodLOC=  MESH.NODES_FACES{isurface} ;  % Nodes pertaining to the surface
DOFrLOC =  small2large(rnodLOC,ndim) ;  % DOFs associatted to the surface. ndim is the number of spatial dimensions
COORface = MESH.COOR(rnodLOC,:) ;   % Spatial coordinates of the nodes 
Centroid = GEOproperties.FACES{isurface}.CENTROID;  % Centroid of the surface
% 4. Coordinate relative to the centroid
COORfaceREL = bsxfun(@minus,COORface',Centroid')' ;  % 



BasisRrb = ConstructBasisRigidBody(COORfaceREL) ; % Rigid body modes 
Rwarping = ConstructWarpingMode(COORfaceREL) ; % Warping mode 
ModesFace = [BasisRrb,Rwarping] ;  

% LOCAL 
% --------  TO BE REVISED, IT DOESN'T WORK
% RotMatrix = [MESH.NormalBoundaryElementsFace{isurface}(:,1),MESH.TangentBoundaryElementsFace{isurface}{1}(:,1),...
%     MESH.TangentBoundaryElementsFace{isurface}{2}(:,1)]  ; 
% COORfaceREL_loc = RotMatrix*COORfaceREL'; 
% COORfaceREL_loc = reshape(COORfaceREL_loc(:),[],size(COORfaceREL,2)) ; 
% BasisRrb_loc = ConstructBasisRigidBody(COORfaceREL_loc) ; % Rigid body modes 
% Rwarping_loc = ConstructWarpingMode(COORfaceREL_loc) ; % Warping mode



nloads = length( LOCVAR.PRESCRIBED_DISP);
% Uloc = zeros(length(DOFrLOC),nloads) ;
% aloc = zeros(nloads,length(DATA.STEPS)) ;
 ntimesteps = length(DATA.STEPS) ;

U = zeros(size(ModesFace,1),nloads) ;
a = zeros(nloads,ntimesteps) ;

for  iload = 1:nloads
     LOC_DISPF  = LOCVAR.PRESCRIBED_DISP(iload) ;
    AMPLITUDE= LOC_DISPF.AMPLITUDE ;  
    INTERVAL = LOC_DISPF.INTERVAL ;    
     TIMEFUN = LOC_DISPF.TIMEFUN ;   
%      LOC_DISPF = DefaultField(LOC_DISPF,'ISLOCAL',0) ; 
%      if LOC_DISPF.ISLOCAL == 1
%          ModesFace = [BasisRrb_loc,Rwarping_loc] ;  
%      end
     
  %  D= zeros(length(DOFrLOC),ntimesteps)   ;
    
  U(:,iload) = ModesFace*AMPLITUDE(:) ;     
  FactorSteps =  TIMEFUN(DATA.STEPS) ;
  FactorSteps = FactorSteps.*(DATA.STEPS >= INTERVAL(1)) ;
  FactorSteps = FactorSteps.*(DATA.STEPS <= INTERVAL(2)) ;
  a(iload,:) = FactorSteps ; 
  
end
% 
% 
% U = cell2mat(U) ;
% a = cell2mat(a) ;




