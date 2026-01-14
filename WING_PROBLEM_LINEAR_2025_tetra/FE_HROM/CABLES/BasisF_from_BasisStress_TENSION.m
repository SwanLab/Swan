function SNAPforceS = BasisF_from_BasisStress_TENSION(BstRED_l,BasisPone,DATA)

if nargin == 0
    load('tmp.mat')
end

% 

nDEF = size(BstRED_l,2) ; % Number of displacement modes
nBasisPone = size(BasisPone,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
% DATA = DefaultField(DATA,'NO_USE_Deformation_gradient_in_Small_Strains',0) ; 
% DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ; 
% if DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1  &&  DATA.SMALL_STRAIN_KINEMATICS ==1  
%     nstrainF = DATA.MESH.nstrain ;  
% else
% nstrainF = DATA.MESH.ndim^2; 
% end
nstrainF = DATA.MESH.ndim ; 

for I = 1:nDEF
    SNAPloc = zeros(DATA.MESH.ngausT,nBasisPone) ;
    for istrain = 1:nstrainF
        B_stress = bsxfun(@times,BasisPone(istrain:nstrainF:end,:),BstRED_l(istrain:nstrainF:end,I));
        % if DATAIN.NOTMULTIPLIED_BY_WEIGHTS ==0
       % SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ;  %
       % Removed 12-APril-2020
        % else
             SNAPloc = SNAPloc + B_stress ;
        % end
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;
end