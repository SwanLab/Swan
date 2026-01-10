function [MESH1Ddom,MESH3D] = GeometryStructure(MESH1Dinput,MESH3D,DATAIN) ;
% See BeamROM.tex, IMPLEM1.tex
%   MESH1D =
%
%                 COOR: [101x3 double]    --> Coordinates of all 1D nodes
%      ELEMENTS_MESH1D: [1x94 double]   --> 1D elements that form the skeleton structure
%          NNODES_elem: [94x1 double] ---> Number of nodes for each 1D entities.  = 2 for slices and two-faces joints. >2.Joints
%                   CN: [94x2 double]  ->> Connectivities
%     NODES_INTERFACES: [95x1 double]   --> Nodes where interfaces are
%     located
%         ROTATIONSdom: [3x282 double]   --> Rotations slices/joints
%              TRANSFM: [1x1 struct]    --> Affine transformation
%         ROTATIONSint: [3x303 double]  --> Rotations interfaces
%                 TYPE: [1x1 struct]   ---> Type of strucutral members
% --------------------------
%%
% MESH1D.TRANSFM  --> MAtrices defining affine transformation of each
% slice/joint
%
% ans =
%
%      A: [3x282 double]
%      D: [3x282 double]
% %     a0: [3x94 double]
% %
% MESH1D.TYPE
%
% ans =
%
%                  ISBEAM: [94x1 double]
%                 INDEX3D: [94x1 double]
%        INDEX_BEAM_JOINT: [94x1 double]
%     SUBINDEX_BEAM_JOINT: [94x1 double]
%----------------------------------------------------------------------------



if exist('ReadMeshFile','file') == 0 ; addpath('FE_CODE') ; end
if exist('GeometryStructure','file') == 0 ; addpath('FE_CODE/BeamROM') ; end
% Typical inputs
% -------------
% Structure array: MESH1D, MESH3D, DATAINs
%

% --------------------------------
% Extracting 3D properties SLICES.
% --------------------------------
disp('----------------------------------------')
disp('  3D properties SLICES')
disp('----------------------------------------')
for islice =1:length(MESH3D.SLICES)
    disp(['islice =',num2str(islice)])
    DATA3D_loc  =  GeometrySlice(MESH3D.SLICES(islice),DATAIN) ;
    
    if ~isempty(DATAIN.SCALE_FACTOR)  
        % Scaling slices in the X-Y and Z direction
        DATA3D_loc = ScalingFactorsSlice(DATA3D_loc,DATAIN.SCALE_FACTOR) ;
    end
    
    DATAIN = DefaultField(DATAIN,'angDOM',[]) ;
    DATAIN = DefaultField(DATAIN,'ELEVATION_Z',[]) ; 
      if ~isempty(DATAIN.ELEVATION_Z)
        % For helicoidal structures
        % -------------------------
        DATA3D_loc = ElevatedDomains(DATA3D_loc,DATAIN) ;
    end
    if ~isempty(DATAIN.angDOM) && DATAIN.angDOM ~= 0
        % Curved elements
        DATA3D_loc = CurvedDomains_beamsGEN(DATA3D_loc,DATAIN.angDOM,DATAIN) ;
    else        
        DATA3D_loc.rotDOMfacesLOC = cell(1,length(DATA3D_loc.CENTRf1)) ;
    end    
  
    MESH3D.SLICES(islice).DATA3D = DATA3D_loc ;    
end
disp('...END')
% --------------------------------
% Extracting 1D properties (except for rotations matrix for joints)
% --------------------------------
disp('----------------------------------------')
disp('  1D properties BEAMS')
disp('----------------------------------------')
[MESH1D,xSIGN ]= Geometry1Dstructure(MESH1Dinput,MESH3D,DATAIN) ;
disp('--------------------------------------------------')
disp('...END')
% --------------------------------
% Extracting 3D properties JOINTS
% % --------------------------------
disp('----------------------------------------')
disp('  3D properties JOINTS')
disp('----------------------------------------')
for itypeJOINT = 1:length(MESH3D.JOINTS)
    disp(['islice =',num2str(itypeJOINT)])
    % ------------------------------------------------------------------------
    % Check first whether the msh file exists (if not, it helps to create it)
    CheckMeshJoint(itypeJOINT,MESH1D,MESH3D,DATAIN ) ;
    % ------------------------------------------------------------------------
    % Reading mesh and geometric entities JOINTS MESH3D.JOINTS(itype)
    LOC =  GeometryJoints(MESH1D,MESH3D,itypeJOINT,DATAIN) ;
    MESH3D.JOINTS(itypeJOINT).DATA3D = LOC ;
end
disp('--------------------------------------------------')
disp('...END')
%----------------------------------------
% Computing rotation matrices for JOINTS. Determining correspondence of
% local numbering and global numbering (in INFOLINES)
%-----------------------------------------
disp('----------------------------------------')
disp('  1D properties JOINTS')
disp('----------------------------------------')
[ROTATIONS,INFOLINES,SUCCESSFUL_MATCH] = LocalAxisDefinitionJoints(MESH1D,MESH3D,DATAIN) ;
MESH1D.ROTATIONS = ROTATIONS;  MESH1D.INFOLINES = INFOLINES ;
MESH1D.SUCCESSFUL_MATCH = SUCCESSFUL_MATCH ;
disp('--------------------------------------------------')
disp('...END')
%-----------------
DATAIN = DefaultField(DATAIN,'PLOT_3D_STRUCTURE_CHECK',0);
% ---------------------
% PRINT MESH 3D DOMAINS
% --------------------
DATAIN.xSIGN = xSIGN ; 
if  DATAIN.PLOT_3D_STRUCTURE_CHECK == 1
    disp('----------------------------------------')
    disp('  PRINTING IN GID')
    disp('----------------------------------------')
    Plot3Dstructure(MESH1D,MESH3D,DATAIN) ;
    disp('--------------------------------------------------')
    disp('...END')
end
%-----------------
%%% Re-writting CONNECTIVITY MATRIX of the 1D mesh, as well as defining
%%% rotation matrices of interfaces
MESH1Ddom = Adapt_MESH1D_DomainDecomposition(MESH1D,MESH3D) ;
MESH1Ddom.SUCCESSFUL_MATCH = SUCCESSFUL_MATCH ;
% -------------------------------------------------------------------------


