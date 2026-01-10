function [RotMatrixGlo,ORDER_CONNECTIVITIES,CN,xSIGN] =RotationMatrixSlices(MESH1D,itypeBEAM,RotMatrixGlo,...
    CENTRf1,CENTRf2,itypeSLICE,ORDER_CONNECTIVITIES,ROT_LOC_SLICE_normals,DATAIN)
% Computation of rotation matrix for 1D mesh (SLICES)
% Routine conceived for: (from 11-Jan-2019)
% 1) Straight beams 
% 2) Curved beams whose directrix lies in the XY plane 
% 3) Beams in which the slice is rotated a
%
% JAHO, 22-January-2018
% ----------------------
if nargin == 0
    load('tmp.mat')
end

% Type of elements associated to beam "itypeBEAM"
INDEXelem = find(MESH1D.MaterialType == itypeBEAM) ;
PROP = MESH1D.PROP(itypeBEAM) ;


% Definition of local z-axis   
% ---------------------------
[ROT_LOC_SLICE,rotX,rr1,e3] =  DefinitionZaxis(PROP,MESH1D,CENTRf1,CENTRf2,DATAIN) ; 




%% Definition of xLOCAL
%-------------------------------------------------------------
[xSIGN] = Xlocal_Slice(PROP,rr1,MESH1D,INDEXelem) ;
% -----------------------------------------------------------------
if xSIGN ==-1
    ORDER_CONNECTIVITIES(:) = -1 ;
end
DATAIN = DefaultField(DATAIN,'ISTWIST_ANGLE',0) ;

if DATAIN.ISTWIST_ANGLE == 1
    % Twisting angle
    rot_TWIST =  ROT_LOC_SLICE_normals ;
    % Minimum X
    xMIN =  min(MESH1D.COOR(:,1)) ;    
else
    rot_TWIST = eye(length(CENTRf1)) ;
    xMIN  = 0 ;    
end



% Loop over elements

for eLOC = 1:length(INDEXelem)
    eGLO = INDEXelem(eLOC) ;
    nodoINI = MESH1D.CN(eGLO,1) ; % Initial node
    nodoFIN = MESH1D.CN(eGLO,2) ; % Final node
    xINI = MESH1D.COOR(nodoINI,:) ; 
    xFIN = MESH1D.COOR(nodoFIN,:) ;
  
    % Rotation matrix
    % First vector
   [rotMATloc,e1] = LocalRotationMatrixNEW(xSIGN,xFIN,xINI,e3,rotX,...
    DATAIN,xMIN,ROT_LOC_SLICE) ; 
    %-----------------------------------------------------------
    ndim = length(e1) ;
    ifin = ndim*eGLO ;
    iini = ifin-ndim+1 ;
    RotMatrixGlo(:,iini:ifin) = rotMATloc;
    
end


CN = MESH1D.CN ;


% %%% NOW WE CHANGE THE ORDER OF THE CONNECTIVITIES
% if any(ORDER_CONNECTIVITIES==-1)
% disp('----------------------------------')
% disp('CHANGING order of connectivities')
% CN = MESH1D.CN ;
%  CN(INDEXelem,:) = [MESH1D.CN(INDEXelem,2),MESH1D.CN(INDEXelem,1)   ] ;
%   MESH1D.PROP(itypeBEAM).xLOCAL = 1 ;  ;
% %%%
%
% % PROP = DefaultField(PROP,'xLOCAL',1) ; %  Direction local x
% % if ~isempty(PROP.xLOCAL)
% %     xSIGN = PROP.xLOCAL ;
% % else
% %     xSIGN = 1;
% % end
%
% end



