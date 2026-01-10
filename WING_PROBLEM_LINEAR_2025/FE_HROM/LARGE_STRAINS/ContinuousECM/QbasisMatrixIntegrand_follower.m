function Q = QbasisMatrixIntegrand_follower(NstRED_l,Basis_tpress,DATA,wSTs,DATAoffline)

if nargin == 0 
    load('tmp.mat')
     
end

SNAPred_FPRESS = BasisF_from_Basis_tpress(NstRED_l,Basis_tpress,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;
SNAPred_FPRESS = bsxfun(@times,SNAPred_FPRESS,sqrt_wST) ;
% Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
%DATAoffline.errorFINT = 1e-3;
DATAsvd.RELATIVE_SVD = 1;
DATAsvd.HIDE_OUTPUT = 1 ; 
[Q,S,V,eSVD,Rsup] = RSVDT(SNAPred_FPRESS,DATAoffline.errorFINT,[],0,DATAsvd) ;

if DATAoffline.errorFINT == 0
    ifig = 3000 ;
    SVDplotERROR_local(S,ifig) ;
end
% % Enlarge the basis matris for SNAPredFINT
a  = sqrt_wST - Q*(Q'*sqrt_wST) ;
if norm(a) > 1e-10
    a = a/norm(a) ;
    Q = [Q,a] ;
end