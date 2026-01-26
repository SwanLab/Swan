function [A,D,a0] = TransformationMatrices(MESH1D,MESH3D,itypeBEAM,itypeSLICESGLO,A,D,a0,jbeam)

if nargin == 0
    load('tmp.mat')
end

% Indices elements and nodes of beam "itypeBEAM", subtype jbeam
INDEXelem =  MESH1D.INFOLINES.ELEMENTS{itypeBEAM}{jbeam}     ; 
INDEXnodes =  MESH1D.INFOLINES.NODES{itypeBEAM}{jbeam}     ;  



PROPS = MESH1D.PROP(itypeBEAM) ; 
PROPS = DefaultField(PROPS,'ELEMENT_OF_EACH_ENTITY',[]) ; 
if isempty(PROPS.ELEMENT_OF_EACH_ENTITY)
     if length(itypeSLICESGLO) == 1
        PROPS.ELEMENT_OF_EACH_ENTITY = {1:length(INDEXelem)} ; 
     else
         error('You must define field "ELEMENT_OF_EACH_ENTITY" ')
    end
end


% GEOMETRY REFERENCE DOMAINS
% --------------------------
GEOref = [] ; 
TypeOfSlices = MESH1D.TypeOfSlices ;  % Type of Slices for beam  itypeBEAM  
for iLOC = 1:length(itypeSLICESGLO)
    itypeSLICE = itypeSLICESGLO(iLOC) ; 
    DATA3D = MESH3D.SLICES(itypeSLICE).DATA3D ;
    
    Cf1 = DATA3D.CENTRf1; %   Centroid face F1, parent slice
    Cf2 = DATA3D.CENTRf2; %   Centroid face F2, parent slice
    %rr1 = (Cf2-Cf1)/norm(Cf2-Cf1) ;  
    %curvature_theta_max  = -asin(rr1(2)) ; % Default correction angle, curved elements
    %CenterRotation = DATA3D.CenterRotation ; 
    
    GEOref(iLOC).dX_0 = norm(Cf2-Cf1) ; % Length slice, reference configuration
    
    %% Coordinates related to centroid Cf1
    COORrel = zeros(size(DATA3D.COOR)) ;
    for idim=1:length(Cf1)
        COORrel(:,idim) = DATA3D.COOR(:,idim) - Cf1(idim) ;
    end
    
    % Maximum y-coordinate , z-coordinate
    BNDnodes = unique(DATA3D.CNb(:))  ;
    COORbnd = COORrel(BNDnodes,:) ;
    ndim = size(COORbnd,2) ;
    GEOref(iLOC).YMAX = max(COORbnd(:,2)) ;
    GEOref(iLOC).YMIN = min(COORbnd(:,2)) ;
    if ndim==3
    GEOref(iLOC).ZMAX = max(COORbnd(:,3)) ;
    GEOref(iLOC).ZMIN = min(COORbnd(:,3)) ;
    end
    
 %  TypeOfSlices(INDEXelem(PROPS.ELEMENT_OF_EACH_ENTITY{iLOC})) = iLOC ; 
    
end

% 1D elements of type "itypeBEAM"
% PROPS = DefaultField(PROPS,'ELEMENT_OF_EACH_ENTITY',[]) ; 


%             MESH1D.PROP(itypeBEAM).CROSS_SECTION.YMAX.FUN = @FSINUSOIDAL   ; % SINUSOIDAL FUNCTION
%             MESH1D.PROP(itypeBEAM).CROSS_SECTION.YMAX.ARGS = [3,0.2]  ; % ARGS(1) = Number of periods, ARGS(2) = Amplitude
%             MESH1D.PROP(itypeBEAM).CROSS_SECTION.YMIN.FUN = @FSINUSOIDAL   ; % SINUSOIDAL FUNCTION
%             MESH1D.PROP(itypeBEAM).CROSS_SECTION.YMIN.ARGS = [3,-0.2]  ; % ARGS(1) = Number of periods, ARGS(2) = Amplitude

% ---------------
% Default values, cross section variations  (constant)
% ---------------
PROPS = DefaultField(PROPS,'CROSS_SECTION',[]) ;
CROSS_SECTION = PROPS.CROSS_SECTION ;

CROSS_SECTION = DefaultField(CROSS_SECTION,'YMAX',[]) ;
CROSS_SECTION.YMAX = DefaultField(CROSS_SECTION.YMAX,'FUN','FCONSTANT') ;
CROSS_SECTION.YMAX = DefaultField(CROSS_SECTION.YMAX,'ARGS',0) ;

CROSS_SECTION = DefaultField(CROSS_SECTION,'YMIN',[]) ;
CROSS_SECTION.YMIN = DefaultField(CROSS_SECTION.YMIN,'FUN','FCONSTANT') ;
CROSS_SECTION.YMIN = DefaultField(CROSS_SECTION.YMIN,'ARGS',0) ;


if ndim == 3
    CROSS_SECTION = DefaultField(CROSS_SECTION,'ZMAX',[]) ;
    CROSS_SECTION.ZMAX = DefaultField(CROSS_SECTION.ZMAX,'FUN','FCONSTANT') ;
    CROSS_SECTION.ZMAX = DefaultField(CROSS_SECTION.ZMAX,'ARGS',0) ;
    
    CROSS_SECTION = DefaultField(CROSS_SECTION,'ZMIN',[]) ;
    CROSS_SECTION.ZMIN = DefaultField(CROSS_SECTION.ZMIN,'FUN','FCONSTANT') ;
    CROSS_SECTION.ZMIN = DefaultField(CROSS_SECTION.ZMIN,'ARGS',0) ;
    
end
 


% Which nodes is the "first one" --- in the sense of the local x-axis? Let us interrogate the matrix of
% rotations
% First of all, we permute the columns of the connectivitie matrix, just in case
PROPS= DefaultField(PROPS,'xLOCAL',1) ; % = 1 ;
SIGNO = PROPS.xLOCAL  ;
firstELEMENT = INDEXelem(1) ;  % First elem
firstNODE = INDEXnodes(1) ;  % First node of this set 
% Which position occupies "firstNODE"
CNfirst = MESH1D.CN(firstELEMENT,:) ; 
if ~isempty(SIGNO) && SIGNO ==-1     
    CNfirst = [CNfirst(2) CNfirst(1)] ;
end
if CNfirst(1) ~= firstNODE
    % we have to invert the order of both lines and nodes 
    NELEM = length(INDEXelem) ; NNODE = length(INDEXnodes)  ; 
    INDEXelem = INDEXelem(NELEM:-1:1) ; 
    INDEXnodes = INDEXnodes(NNODE:-1:1) ; 
end

%%% Now we have to compute the curvilinear length of the beam 
% -----------------------------------------------------------


 
sCURV = zeros(length(INDEXelem),2) ; % Curvilinear coordinates
sINI = 0;  
% ----------------------------------------------------
for eLOC = 1:length(INDEXelem)
    nodoINI = INDEXnodes(eLOC) ;
    nodoFIN = INDEXnodes(eLOC+1) ;
    xINI = MESH1D.COOR(nodoINI,:) ;
    xFIN = MESH1D.COOR(nodoFIN,:) ;
    % Curvilinear coordinate
    sFIN = sINI +norm(xFIN-xINI) ;
    sCURV(eLOC,:) = [sINI,sFIN] ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sINI = sFIN ;  % Update of curvilinear coordinate
end


% Computing matrices A,D and a0, see BeamROM.pdf
LENGTH_STRUCTURE =  sCURV(end,2) ;

for eLOC = 1:length(INDEXelem)
    eGLO = INDEXelem(eLOC) ;
    if length(itypeSLICESGLO) == 1
        iSLICE = 1; 
    else
    iSLICE = TypeOfSlices(eGLO) ; 
    iSLICE = find(itypeSLICESGLO==iSLICE) ; 
    end
    
    [A_loc,D_loc,a0_loc] = MatricesVaryingCrossSection(GEOref(iSLICE),CROSS_SECTION,sCURV(eLOC,:),...
        LENGTH_STRUCTURE) ;
    a0(:,eGLO) = a0_loc;
    ifin = eGLO*ndim ;
    iini = eGLO*ndim-ndim+1 ;
    A(:,iini:ifin) = A_loc;
    D(:,iini:ifin) = D_loc;   
end
