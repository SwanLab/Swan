function [SNAPforceS] = BasisFfromStress_solo(BasisS,BdomRED, Wdom)
% See Implementation.pdf, section MAtrix of internal forces
%dbstop('4')
if nargin ==0
    load('tmp1.mat')
   % DATACUB = [] ;
end

ngaus = length(Wdom);
nstrain = size(BdomRED,1)/length(Wdom); 
nDEF = size(BdomRED,2) ; % Number of displacement modes
nBASISs = size(BasisS,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
for I = 1:nDEF
    SNAPloc = zeros(ngaus,nBASISs) ;
    for istrain = 1:nstrain
        B_stress = bsxfun(@times,BasisS(istrain:nstrain:end,:),BdomRED(istrain:nstrain:end,I)); 
        SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ; 
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;     
end

% 
% [U,S,V,h3,h4] = SVD_and_error(SNAPforceS,nfigure,LEGENDG,[],COLOR,DATAIN ) ;
% 


% 
% if ~exist('RSVDT')
%     addpath('SVDlibrary')
% end
% 
% 
% 
% 
% hold on
% nfigure = 100;
% LEGENDG = 'Internal Virtual Work ' ;
% COLOR = 'k-' ;
% %dbstop('40')
% DATAIN = DefaultField(DATAIN,'TOL_LOC_InternalForces',1e-6) ; 
% DATAIN.TOL_LOC = DATAIN.TOL_LOC_InternalForces  ; 
% [U,S,V,h3,h4] = SVD_and_error(SNAPforceS,nfigure,LEGENDG,[],COLOR,DATAIN ) ;
% hold on
% legend([h3 ],{'Int. virt. work'})
% legend([h4 ],{'Int. virt. work'})
% 
% BasisF = U  ; % = [2] ;
% SingVal_intf = S ;  
% % save(DATA.NAMEWS,'-append','BasisS','SingVal_intf');
% 
% 
