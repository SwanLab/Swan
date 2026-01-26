function [Q,SNAPredFINT_nw] = QbasisMatrixIntegrand(BstRED_l,BasisPone,DATA,wSTs,DATAoffline)

if nargin == 0
    load('tmp.mat')
    % DATAoffline.errorFINT = 1e-6 ;
end

SNAPredFINT_nw = BasisF_from_BasisStress_PK1(BstRED_l,BasisPone,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;
SNAPredFINT = bsxfun(@times,SNAPredFINT_nw,sqrt_wST) ;
% Determine an  orthogonal basis matrix $Q$ for the column space of $\SNAPredFINTw{}{}$
%DATAoffline.errorFINT = 1e-3;


DATAoffline = DefaultField(DATAoffline,'USE_SVD_RANDOMIZED_FOR_FINT',1) ; % = 0 ;
DATAsvd.RELATIVE_SVD = 1;
DATAsvd.HIDE_OUTPUT = 1 ;
if DATAoffline.USE_SVD_RANDOMIZED_FOR_FINT == 1
    
    if norm(SNAPredFINT,'fro') == 0
        Q = [] ; 
    else
    [Q,S,V,eSVD,Rsup] = RSVDT(SNAPredFINT,DATAoffline.errorFINT,[],0,DATAsvd) ;
    end
    
    
    
else
    [Q,S,V,eSVD] = SVDT(SNAPredFINT,DATAoffline.errorFINT,DATAsvd) ;
    
end
if DATAoffline.errorFINT == 0
    ifig = 3000 ;
    plot_svd_error = 0 ;
    if plot_svd_error == 1
        SVDplotERROR_local(S,ifig) ;
    end
end

% % Enlarge the basis matris for SNAPredFINT

if isempty(Q)
    Q =  sqrt_wST/norm(sqrt_wST) ;
else
    a  = sqrt_wST - Q*(Q'*sqrt_wST) ;
    if norm(a) > 1e-10
        a = a/norm(a) ;
        Q = [Q,a] ;
        SNAPredFINT_nw = [ones(size(a)),SNAPredFINT_nw] ; %
        
    end
end