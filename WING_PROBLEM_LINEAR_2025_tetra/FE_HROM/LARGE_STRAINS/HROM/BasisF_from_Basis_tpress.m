function SNAPforceS = BasisF_from_Basis_tpress(NstRED_l,Basis_tpress,DATA)

if nargin == 0
    load('tmp1.mat')
end

% 

nDEF = size(NstRED_l,2) ; % Number of displacement modes
nBasist_press = size(Basis_tpress,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
nstrainF = DATA.MESH.ndim; 

for I = 1:nDEF
    SNAPloc = zeros(DATA.MESH.HYDRO.ngausT,nBasist_press) ;
    for istrain = 1:nstrainF
        B_stress = bsxfun(@times,Basis_tpress(istrain:nstrainF:end,:),NstRED_l(istrain:nstrainF:end,I));
        % if DATAIN.NOTMULTIPLIED_BY_WEIGHTS ==0
       % SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ;  %
       % Removed 12-APril-2020
        % else
             SNAPloc = SNAPloc + B_stress ;
        % end
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;
end