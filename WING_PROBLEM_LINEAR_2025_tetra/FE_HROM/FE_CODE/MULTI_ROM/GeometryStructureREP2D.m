function [MESH2Ddom,MESH3D] = GeometryStructureREP2D(MESH2Dinput,MESH3D,DATAIN) ;
% Copy of GeometryStructure.m 
% See BeamROM.tex, IMPLEM1.tex
%   MESH2D =
%
%                 COOR: [101x3 double]    --> Coordinates of all 1D nodes
%      ELEMENTS_MESH2D: [1x94 double]   --> 1D elements that form the skeleton structure
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
% MESH2D.TRANSFM  --> MAtrices defining affine transformation of each
% slice/joint
%
% ans =
%
%      A: [3x282 double]
%      D: [3x282 double]
% %     a0: [3x94 double]
% %
% MESH2D.TYPE
%
% ans =
%
%                  ISBEAM: [94x1 double]
%                 INDEX3D: [94x1 double]
%        INDEX_BEAM_JOINT: [94x1 double]
%     SUBINDEX_BEAM_JOINT: [94x1 double]
%----------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end


% Typical inputs
% -------------
% Structure array: MESH2D, MESH3D, DATAINs
%

% --------------------------------
% Extracting 3D properties SLICES.
% --------------------------------
disp('----------------------------------------')
disp('  3D properties RVES')
disp('----------------------------------------')
DATAIN = DefaultField(DATAIN,'angDOM',[]) ; % Angle formed by faces f1 and f3
for irve =1:length(MESH3D.RVES)
    disp(['irve =',num2str(irve)])
    DATA3D_loc  =  GeometryRVE(MESH3D.RVES(irve),DATAIN) ;
    if ~isempty(DATAIN.angDOM) && DATAIN.angDOM ~= 0
        % Curved elements
        MESH3D.RVES(irve).DATA3D  = CurvedDomains(DATA3D_loc,DATAIN.angDOM) ;
    else
        % No curved elements 
        MESH3D.RVES(irve).DATA3D  = DATA3D_loc ;        
        MESH3D.RVES(irve).DATA3D.rotDOMfacesLOC = cell(1,length(DATA3D_loc.COOR_MIDSIDE_FACE)) ;
    end
    
end
disp('...END')
% --------------------------------
% Extracting 1D properties (except for rotations matrix for joints)
% --------------------------------
disp('----------------------------------------')
disp('  2D properties ')
disp('----------------------------------------')
MESH2D = Geometry2Dstructure(MESH2Dinput,MESH3D,DATAIN) ;
disp('--------------------------------------------------')
disp('...END')


DATAIN = DefaultField(DATAIN,'PLOT_3D_STRUCTURE_CHECK',1);

% ---------------------
% PRINT MESH 3D DOMAINS
% --------------------
if  DATAIN.PLOT_3D_STRUCTURE_CHECK == 1
    disp('----------------------------------------')
    disp('  PRINTING IN GID')
    disp('----------------------------------------')
    Plot3DstructureRVE(MESH2D,MESH3D,DATAIN) ;
    disp('--------------------------------------------------')
    disp('...END')
end
 
% -------------------------------------------------------------------------
MESH2Ddom = MESH2D ;

for ientity = 1:length(MESH2Ddom.PROP)
    switch  MESH2Ddom.PROP(ientity).TYPE
        case 'RVE'
            islice= MESH2Ddom.PROP(ientity).INDEX_ENTITY ;
            MESH2Ddom.PROP(ientity).NameWSmodes = MESH3D.RVES(islice).NameWSmodes ;
        otherwise
            error('Option not implemented')
    end
end

