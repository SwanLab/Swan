function [Vrb,Vdef] = CandidatesInterfaceModesRVE(BasisUdef,f1,f2,M,Vrb,DATAIN)
%  Copy of CandidatesInterfaceModes (this was made for beams, while the present rouitine is for RVEs)
if nargin ==0
    load('tmp2.mat')
end

BasisREFERENCE = [BasisUdef(f1,:), BasisUdef(f2,:)];
%%% PURGE ROTATION MODES 
% ----------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix Vrb
Mchol = chol(M) ;
nmodesINTF = size(BasisUdef,2) ;
Vdef = MorthogonalMatrix(M,Mchol,BasisREFERENCE,Vrb,DATAIN)  ;
% if size(Vdef,2) > nmodesINTF 
% BasisREFERENCE = BasisREFERENCE(:,1:nmodesINTF) ; % There cannot be more interface modes than disp. modes
% end
%%%%%%%%%%%%%%%
% Now we turn Vrb M-orthogonal
Xbar = Mchol*Vrb ; 
TOL = 1e-6 ; 
[Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
Vrb = Mchol\Ubar ;

 