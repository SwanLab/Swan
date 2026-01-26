function [INFO_RVE,MESHdom,OPERFE,MATPRO]= IndexesREFdomain(IDOM,IndexBoundaries,MESH,DATA,OPERFE,MATPRO)
%--------------------------------------------------------------------------
% FUNCTION: IndexesREFdomain
%
% PURPOSE:
%   Extracts the indices and finite element structures corresponding to a
%   specific subdomain within a global mesh. It returns the reduced mesh
%   (`MESHdom`), finite element operators (`OPERFE`), and material properties
%   (`MATPRO`) restricted to the subdomain labeled by `IDOM`.
%
%   This is a key step in multiscale and localized reduced-order models (e.g.,
%   EIFEM), where operators and snapshots must be localized to a region of
%   interest (ROI), such as a representative volume element (RVE).
%
% INPUTS:
%   - IDOM            : Identifier of the target subdomain (integer index)
%   - IndexBoundaries : Array of boundary label indices to extract boundary nodes
%   - MESH            : Global mesh structure with fields like COOR, CN, etc.
%   - DATA            : Data structure containing mesh metadata (e.g., ndim, ngaus)
%   - OPERFE          : Global finite element operator structure (may be empty)
%   - MATPRO          : Global material property matrices (e.g., constitutive tensors)
%
% OUTPUTS:
%   - INFO_RVE        : Structure with index maps (nodes, DOFs, elements, Gauss pts)
%   - MESHdom         : Submesh structure corresponding to the domain `IDOM`
%   - OPERFE          : Restricted FE operator structure for the subdomain
%   - MATPRO          : Restricted material property fields for the subdomain
%
% FUNCTIONAL OVERVIEW:
%   1. Identifies the global nodes and DOFs of the subdomain using `MESH.NODES_DOM{IDOM}`.
%   2. Determines which elements are associated with those nodes (BODYELEM_globNUM).
%   3. Computes the Gauss point indices associated with those elements, for:
%        - scalar variables (e.g., energy densities),
%        - vector/tensor variables (e.g., stress components).
%   4. Constructs the reduced mesh `MESHdom`:
%        - Coordinates,
%        - Connectivity (renumbered locally),
%        - Material types and element types,
%        - Boundary connectivities (if `IndexBoundaries` is given).
%   5. Extracts subdomain-specific FE operators from `OPERFE`:
%        - Strain-displacement matrix B (`Bst`),
%        - Integration weights (`wSTs`).
%   6. Truncates material properties in `MATPRO` to the selected Gauss points or elements,
%      depending on the size of each field.
%
% USAGE CONTEXT:
%   Called during the offline stage of multiscale ROM or hyperreduction setups
%   to isolate a subdomain for projection, basis construction, or integration point
%   selection.
%
% NOTES:
%   - Supports both internal and boundary elements.
%   - Mapping and renumbering ensure compatibility with local basis construction.
%   - Compatible with small-strain formulations (extensions to large strain require updates).
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 1-Feb-2023
%   Modified: 23-Sep-2023 to include support for boundary face mapping
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
end

% This function returns INFO_RVE.{PROPERTY} for the selected subdomain
% IDOM (PROPERTY may be the indexes of the  nodes,  the dOFS, GAUss points indexes, the mesh...etc)
% See for instance: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/01_INNER_DOM/OFFLINE_STEPSfun.m
% JAHO 1-Feb-2023
% ----------------------------------------------
% NODES OF the DOMAIN
% --------------------------------------------------------
INFO_RVE.NODES_globNUM = MESH.NODES_DOM{IDOM} ;  %  The are the numbering of the nodes of the selected domain
% DOFs
ndim = DATA.MESH.ndim ; % Number of spatial dimensions
INFO_RVE.DOFS_globNUM = small2large(INFO_RVE.NODES_globNUM,ndim); % DOFs associated to such nodes
% ELEMENTS ASSOCIATED WITH THESE NODES
% ---------------------------------------------------------
[dummy,  INFO_RVE.BODYELEM_globNUM]= ElemBnd(MESH.CN,INFO_RVE.NODES_globNUM); % Here we identify the elements of the domain
% ---------------------------------------------------------
ngaus_STRESS = DATA.MESH.ngaus_STRESS ; % Number of Gauss points integration of stresses
INFO_RVE.GaussINDEX_scalarSTR = small2large(INFO_RVE.BODYELEM_globNUM,ngaus_STRESS);  % Indexes Gauss points for scalar variables (Left-Hand side)
nstrain = DATA.MESH.nstrain; % Number of strain/stress variables
INFO_RVE.GaussINDEX_stress = small2large(INFO_RVE.GaussINDEX_scalarSTR ,nstrain);

% FINITE ELEMENT/GEOMETRIC  OPERATORS OF THE SELECTED MESH **+
% ------------------------------------------------------------
% COORDINATES
MESHdom.COOR = MESH.COOR(INFO_RVE.NODES_globNUM,:) ;

if isfield(MESH,'INFO_PERIODIC_CONDITIONS')
  MESHdom.INFO_PERIODIC_CONDITIONS = MESH.INFO_PERIODIC_CONDITIONS ; 
end

% ---------------------------------------------------
% BODY ELEMENTS
% CONNECTIVITIES
CN_GLO = MESH.CN(INFO_RVE.BODYELEM_globNUM,:) ;  % Table of connect. Global numbering 
% REnumbering
NODESnew = 1:length(INFO_RVE.NODES_globNUM) ;  
MESHdom.CN = RenumberConnectivities(CN_GLO,NODESnew) ; %Renumbering to LOCAL 
% Material type
MESHdom.MaterialType  = MESH.MaterialType(INFO_RVE.BODYELEM_globNUM) ;
% Type of element
MESHdom.TypeElement = MESH.TypeElement ;


% ----------------------------------
% BOUNDARY  ELEMENTS
% -------------------------------------------

% --------------------------------------------------------------


% NEW BOUNDARY IS FORMED
% ------------------------------------------------------------------------------------------------------
if isempty(IndexBoundaries)
    MESHdom.NODES_FACES = [] ;
    MESHdom.CNb  = []  ;
    MESHdom.TypeElementB  = [] ;
else
    NODES_FACES_glo = MESH.NODES_FACES(IndexBoundaries); % These are the nodes of the interface boundaries 
    MESHdom.NODES_FACES = cell(size(NODES_FACES_glo)) ;
    % Renumbering
    CNb_dom_glo = cell(size(NODES_FACES_glo)) ;
    for ifaces = 1:length(NODES_FACES_glo)
        [~,MESHdom.NODES_FACES{ifaces},~] =     intersect(INFO_RVE.NODES_globNUM,NODES_FACES_glo{ifaces}) ;
        [CNb_dom_glo_iface,dummy  ]= ElemBnd(MESH.CNb,NODES_FACES_glo{ifaces});
        % Now we renumber CNb_dom_glo_iface
        CNb_dom_glo{ifaces} = RenumberConnectivities(CNb_dom_glo_iface,MESHdom.NODES_FACES{ifaces}) ;
    end
    MESHdom.CNb = CNb_dom_glo ;
    MESHdom.TypeElementB = MESH.TypeElementB ;
end
% ------------------------------------------------------------------------------------------------------------


% OPERATORS   (small strains so far, 8-Feb-2023)
% ----------
if ~isempty(OPERFE)
    

    OPERFE.Bst = OPERFE.Bst(INFO_RVE.GaussINDEX_stress,INFO_RVE.DOFS_globNUM) ;
    OPERFE.M = [] ;
    OPERFE.wSTs = OPERFE.wSTs(INFO_RVE.GaussINDEX_scalarSTR) ;
    
    % ---------------------
    % MATERIAL PROPERTIES
    % ---------------------
    fff = fieldnames(MATPRO);
    
    for i = 1:length(fff)
        NameProp = fff{i}  ;
        LENGTH_prop = size(MATPRO.(NameProp),1) ;
        if  LENGTH_prop ==  (DATA.MESH.ngausT)*DATA.MESH.nstrain
            %         if isempty(ECMdata.setPoints)
            %             % MATPRO.(NameProp) = MATPRO.(NameProp)(ECMdata.setElements,:) ;
            %             % % CECM, ERROR FOUND  in 10-DEC-2021
            %
            %             MATPRO.(NameProp) = MATPRO.(NameProp)(setPointsElement,:) ;  % CECM, change 10-Dec-2021
            %         else
            MATPRO.(NameProp) = MATPRO.(NameProp)(INFO_RVE.GaussINDEX_stress,:) ;  % DECM
            %   end
        elseif LENGTH_prop == (DATA.MESH.ngausT)
            MATPRO.(NameProp) = MATPRO.(NameProp)(INFO_RVE.GaussINDEX_scalarSTR,:) ;
        else
            MATPRO.(NameProp) = MATPRO.(NameProp)(INFO_RVE.BODYELEM_globNUM,:) ;
        end
    end
    
end
