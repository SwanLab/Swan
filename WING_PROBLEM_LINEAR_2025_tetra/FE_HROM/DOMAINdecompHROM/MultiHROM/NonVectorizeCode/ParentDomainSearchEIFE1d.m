function [CNnew,TRANSF_COORD,Vall_rot] = ParentDomainSearchEIFE1d(EIFEoper_all,COOR,CNold,TypeElement,DATA)
% -------------------------------------------------------------------------
% Adaptation (simplification) of ParentDomainSearchEIFE to 1D problems 
% JAHO, 13-May-2204
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
% JAHO, 10-March-2023/11-March-2023
% -------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
% ---------------------------------------------------------------------
% DETERMINING PARENT DOMAIN FROM A SET OF POTENTIAL CANDIDATES
% ---------------------------------------------------------------------
 % When data is provided this way, we have to check which transformation
% yields zero dilatation
DATA = DefaultField(DATA,'ToleranceDeformationalPartTransformationParentDomain',1e-5) ;


CNnew = CNold ;
Xe = COOR(CNnew,:)' ;
[TRANSF_COORD] = PolarDecompEIFEperm1D(Xe,EIFEoper_all,DATA) ;
Vall_rot = EIFEoper_all.MODES.Vall ;


  