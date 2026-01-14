function SNAPforceS = BasisRED_BasisStress(BdomRED,BasisS,nstrain,ngaus)

nDEF = size(BdomRED,2) ; % Number of displacement modes
nBASISs = size(BasisS,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
for I = 1:nDEF
    SNAPloc = zeros(ngaus,nBASISs) ;
    for istrain = 1:nstrain
        B_stress = bsxfun(@times,BasisS(istrain:nstrain:end,:),BdomRED(istrain:nstrain:end,I));
        % if DATAIN.NOTMULTIPLIED_BY_WEIGHTS ==0
       % SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ;  %
       % Removed 12-APril-2020
        % else
             SNAPloc = SNAPloc + B_stress ;
        % end
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;
end