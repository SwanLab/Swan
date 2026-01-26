function [CNnew,TRANSF_COORD,PERMUT_chosen] = ParentDomainSearchEIFEcorot(EIFEoper_all,COOR,CNold,TypeElement,DATA)
%--------------------------------------------------------------------------
%  ParentDomainSearchEIFEcorot
%
%  This function identifies the most suitable parent domain (from a set of
%  trained EIFEM candidates) for a given coarse-scale element undergoing
%  large rotations and small strains, using a **corotational formulation**.
%  It searches for the transformation (rotation, translation, scaling) 
%  that minimizes the deformational (dilatational) part of the coordinate 
%  mapping between the physical element and the reference domain.
%
%  The function is an adaptation of `ParentDomainSearchEIFE.m`, specifically
%  tailored to handle rigid-body motions and preserve geometric consistency
%  for high-fidelity variational formulations, as required in the co-rotational
%  EIFEM setting (see Sections 6.3 and 12.8 of *EIFEM_largeROTfinal.pdf*).
%
%  INPUTS:
%    - EIFEoper_all : array of trained EIFEM operators (parent domain candidates)
%    - COOR         : nodal coordinates
%    - CNold        : original connectivity of the element
%    - TypeElement  : type of finite element (e.g., 'Quadrilateral', 'Hexahedra')
%    - DATA         : structure with user options and tolerances
%
%  OUTPUTS:
%    - CNnew        : reordered node indices of the element for correct orientation
%    - TRANSF_COORD : structure containing:
%        > ROTATION: Rotation matrix Q_ini
%        > TRANSLATION: Translation vector to parent domain
%        > SCALEFACTOR: Element-wise scaling parameter
%        > detJe: Jacobian determinant of the transformation
%        > IndexParentDomain: index of the selected candidate
%    - PERMUT_chosen: selected permutation from the list DATA.PERMUT
%
%  The function works by:
%    1. Iterating over candidate parent domains.
%    2. Applying all allowed permutations of the element’s node ordering.
%    3. Computing the polar decomposition and evaluating the norm of the 
%       deformational part of the transformation.
%    4. Selecting the permutation and parent domain with the minimum norm,
%       under a user-specified threshold.
%
%  The selected transformation ensures:
%    - Positive Jacobian determinant
%    - Minimal dilatation error
%    - Correct alignment of reference modes (for later use in downscaling)
%
%  The computed `TRANSF_COORD` is later used to rotate interface modes and 
%  construct geometric and variational operators in the corotational regime.
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Versions: 10-Mar-2023, 11-Mar-2023, 22-Oct-2024
%    Comments by ChatGPT4, 13-May-2025
%
%  Related functions:
%    - PolarDecompEIFEpermCOROT.m
%    - CriterionChooseParentCOROT.m
%    - B_N_matricesEIFEbubCOROT_LRss.m
%
%--------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Copy of ParentDomainSearchEIFE.m, adapted to Co-rotational formulation 
% 
%  Given a EIF element formed by nodes CNold, and the properties of
%  a set of parent domain candidates EIFEoper_all, this function searches
%  for the parent domain that minimizes the dilatational component of the
%  transformation of coordinates from the physical domain to the parent
%  domain
% OUTPUT:
% CNnew: New order of the nodes of the EIF element (permutation)
% TRANSF_COORD =
%
%   struct with fields:
%
%        ROTATION: [2×2 double]
%     TRANSLATION: [2×1 double]
%     SCALEFACTOR: 0.0500
%           detJe: 0.0025
%
% Vall_rot --> Rotated interface displacement modes
% ----------------------------------------------------------
%
% JAHO, 10-March-2023/11-March-2023/22-Oct-2024
% -------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
% ---------------------------------------------------------------------
% DETERMINING PARENT DOMAIN FROM A SET OF POTENTIAL CANDIDATES
% ---------------------------------------------------------------------
icand = 1;  % EIFEoper_all(1), EIFEoper_all(2) ....
% When data is provided this way, we have to check which transformation
% yields zero dilatation
DATA = DefaultField(DATA,'ToleranceDeformationalPartTransformationParentDomain',1e-5) ;
while icand <=length(EIFEoper_all)
    
    
    % PERMUTATION CONNECTIVITIES
    
    TRANSF_COORD_perm = cell(1,length(DATA.PERMUT));
    normDmin = 1e20 ; 
    for   iPERM = 1:length(DATA.PERMUT)
        %   disp(['PERM = ',num2str(iPERM)])
        CNnew = CNold(DATA.PERMUT{iPERM}) ;
       % CNnew_input = CNnew(1:DATA.nnodeE_geometry) ; % 26-Apr-2204
        Xe = COOR(CNnew,:)' ;
        [TRANSF_COORD_perm{iPERM},normD] = PolarDecompEIFEpermCOROT(Xe,EIFEoper_all(icand),DATA) ;
        normDmin = min(normDmin,normD) ; 
    end
    
    [TRANSF_COORD,PERMUT_chosen,CNnew] = CriterionChooseParentCOROT(TRANSF_COORD_perm,CNold,DATA,DATA.PERMUT) ;
    
    
  
    
    if  isempty(TRANSF_COORD)
        icand = icand + 1;
    else
        %EIFEoper = EIFEoper_all(icand);
        TRANSF_COORD.IndexParentDomain = icand ;
       % [ndim,~] = size(Xe) ;
       % Vall_rot = zeros(size( EIFEoper_all(icand).MODES.Vall)) ;
        
        %  UNIFORM_SCALING_REFERENCE_ELEMENT= 1;
        
        %         if  DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
        %             FactorV = (1/TRANSF_COORD.SCALEFACTOR) ;
        %             % THIS IS TO REMIND US THAT THIS SCALING FACTOR ACTUALLY
        %             % APPEARS IN DERIVING THE FE B-matrix. We incorporate it here
        %             % because this version (11-March-2023) can only accurately handle uniform
        %             % dilations
        %         else
        %             FactorV = 1;
        %         end
        
%         for imodes = 1:size( EIFEoper_all(icand).MODES.Vall,2)
%             % Rotation of the interface modes
%             % LOCV = FactorV*(TRANSF_COORD.ROTATION*reshape(EIFEoper_all(icand).MODES.Vall(:,imodes),ndim,[])) ;
%             % disp('Temporal...Erase it')
%             
%             
%             LOCV = (TRANSF_COORD.ROTATION'*reshape(EIFEoper_all(icand).MODES.Vall(:,imodes),ndim,[])) ;
%             Vall_rot(:,imodes) = LOCV(:) ;
%         end
        
        
        
        break
    end
    
    
    
    
end
if icand >length(EIFEoper_all)
  
    error(['THIS ELEMENT HAS NO PARENT DOMAIN; maximum gap =  ',num2str(normDmin),  '(you may fix it by setting  DATAOUT.ToleranceDeformationalPartTransformationParentDomain above this value   )'])
end