function [VAR_SMOOTH_FE,DATA,DATA_ECM] = DataPreparationCECM(MESH,DATA,DATA_ECM,wFE,xFE,VSinv)
%--------------------------------------------------------------------------
% function [VAR_SMOOTH_FE, DATA, DATA_ECM] = DataPreparationCECM(MESH, DATA, DATA_ECM, wFE, xFE, VSinv)
%
% PURPOSE:
%   Prepares the finite element (FE) data structures and auxiliary variables needed 
%   for the Continuous Empirical Cubature Method (CECM). This includes interpolation 
%   order definition, element type handling, and mesh geometry storage. It also sets 
%   up the required input structures for subsequent local polynomial interpolations 
%   and weight optimization stages in CECM.
%
% INPUTS:
%   - MESH     : Structure containing FE mesh information:
%                 * COOR         : Node coordinates
%                 * CN           : Element connectivities
%                 * TypeElement  : Element type ('Quadrilateral', 'Triangle', etc.)
%                 * ngausE       : Number of Gauss points per element
%   - DATA     : General configuration structure. Fields optionally defined:
%                 * posgp_given  : (optional) Positions of Gauss points (unused here)
%                 * xLIM         : Domain bounding box (for visualization)
%                 * Integrand    : Data structure with integrand metadata
%   - DATA_ECM : ECM/CECM-specific structure to be updated
%   - wFE      : Vector of Gauss weights (size = total number of Gauss points)
%   - xFE      : Coordinates of Gauss points (not used explicitly here)
%   - VSinv    : Projection matrix (transpose of right singular vectors / singular values)
%
% OUTPUTS:
%   - VAR_SMOOTH_FE : Structure containing data used during the CECM step:
%                      * ngausE                : Number of Gauss points per element
%                      * COOR, CN              : Mesh coordinates and connectivities
%                      * TypeElement           : Type of finite elements used
%                      * DELTRIANG             : Delaunay triangulation (for plotting/diagnostics)
%                      * IND_POLYG_ELEMENT     : Corner node ordering for each element type
%                      * ORDER_POLYNOMIALS     : Degree of polynomial interpolation per dimension
%                      * wSTs                  : Gauss weights (wFE)
%                      * VSinv                 : Inverse basis matrix for projection in CECM
%
%   - DATA           : Updated with default fields: 'xLIM', 'Integrand'
%   - DATA_ECM       : Updated with Integrand and xLIM fields, passed through from DATA
%
% METHOD:
%   1. Transfer mesh info (coordinates, connectivities, element types).
%   2. Infer polynomial interpolation order from number of Gauss points and mesh dimensionality.
%   3. Define default fields in DATA and DATA_ECM if not already present.
%   4. (If 2D or 3D) Define Delaunay triangulation and polygon node ordering 
%      for later use in visualization or neighbor search.
%
% REMARKS:
%   - This function sets up the "VAR_SMOOTH_FE" structure required for the 
%     interpolatory polynomial reconstruction in the CECM weight-reduction phase.
%   - The Delaunay triangulation is optional but useful for visualizations and local searches.
%   - Polygon node ordering is element-type dependent and assumed to follow Gmsh-like conventions.
%
% REFERENCES:
%   - Hernández et al. (2024), "CECM: A continuous empirical cubature method..." 
%     *Comput. Methods Appl. Mech. Engrg.*, Vol. 418, 116552.
%
% SEE ALSO:
%   - ContinuousECM2023
%   - GenerateLocalInterpolants
%   - VariousPLOTS_CECM2023
%
% AUTHOR:
%   J.A. Hernández (UPC/CIMNE), preprocessing stage for continuous cubature method (2024)
%--------------------------------------------------------------------------


% Printing reduced set of elements
%--------------------------------------
%DATA_DECM = ECMpointsPRINT(DATA_DECM,MESH,HYPERREDUCED_VARIABLES) ;
% -------------------------------------
VAR_SMOOTH_FE.ngausE = MESH.ngausE;
% Matrix of coordinates/connectivities
% ---------------------------
VAR_SMOOTH_FE.COOR = MESH.COOR ;
VAR_SMOOTH_FE.CN = MESH.CN ;
% Type of element
% --------------------------
VAR_SMOOTH_FE.TypeElement = MESH.TypeElement ;
% Dela. triangulation
% --------------------------
if size(MESH.COOR,2) >1
    VAR_SMOOTH_FE.DELTRIANG = delaunayTriangulation(MESH.COOR);
    switch VAR_SMOOTH_FE.TypeElement
        case 'Quadrilateral'
            IND_POLYG_ELEMENT = [1 2 3 4 1] ;
            %   ORDER_POLYNOMIALS  =[norder_poly,norder_poly] ;
        case 'Hexahedra'
            IND_POLYG_ELEMENT = [1 2 3 4 5 6 7 8 1] ;
        case 'Triangle'
            IND_POLYG_ELEMENT = [1 2 3 1] ;
        otherwise
            error('element not implemented')
    end
    VAR_SMOOTH_FE.IND_POLYG_ELEMENT = IND_POLYG_ELEMENT;   % Local numbering of corner nodes (polygon)    
end
% ORDER POLYNOMIALS
% ----------------------------
ndim = size(MESH.COOR,2) ;
DATA = DefaultField(DATA,'posgp_given',[]) ;
norder_poly = round(MESH.ngausE^(1/ndim)-1);
ORDER_POLYNOMIALS = norder_poly*ones(1,ndim);

VAR_SMOOTH_FE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
VAR_SMOOTH_FE.wSTs =  wFE ;
VAR_SMOOTH_FE.VSinv = VSinv ; 

% DATA_ECM = DefaultField(DATA_ECM,'PLOT_INTERNAL_FORCE_MODES',0) ;
% 
% if DATA_ECM.PLOT_INTERNAL_FORCE_MODES == 1
%     DATALOCfint.NameFileMesh = DATA_ECM.NameFileMesh_FINT ;
%     DATALOCfint.MaterialType= MESH.MaterialType ;
%     GidPostProcessModesFINT(MESH.COOR,MESH.CN,...
%         MESH.TypeElement,HYPERREDUCED_VARIABLES.PHI,DATA.MESH.posgp,DATALOCfint);
% end

DATA= DefaultField(DATA,'xLIM',[]) ; 
DATA = DefaultField(DATA,'Integrand',[]) ; 
DATA_ECM.Integrand = DATA.Integrand ; 
DATA_ECM.xLIM = DATA.xLIM ; 