function [qDEF,qRB,rDEF,rRB,pINT] = EquationsAssemblyInterfaceM(Tcond,Tr,cBC_d,cBC_r,KdomREDglo,...
    Hdr,Hrd,Hrr,Hcond,fextRED_d,fextRED_r) ; 

if nargin == 0
    load('tmp2.mat') 
end


%%%%%%%%%%%%%%%%%%
% LEFT-HAND SIDE %
%%%%%%%%%%%%%%%%%%

% 1) Matrix Arr (flexibility)
% ---------------------------
% \Amat{rr} =  {\Hcond{T}}\KdomRED{-1}  \Hcond{} 
Arr = Hcond'*(KdomREDglo\Hcond) ;

% 2) Interface stiffness matrix 
% -----------------------------
% \Kintf \defeq \Amat{pr}\Amat{rr}^{-1}\Amat{rp} = \Tcond{} ({\Hcond{T}}\KdomRED{-1}  \Hcond{})^{-1} \Tcond{T}
Kintf = Tcond*(Arr\Tcond') ; 

%%%%%%%%%%%%%%%%%%
% RIGHT-HAND SIDE %
%%%%%%%%%%%%%%%%%%

% 3) Condensed external forces 
% \fextREDcond \defeq \fextRED{}{}_d - \Hqr{}_{dr} \Hqr{-1}_{rr} \fextRED{}_r
fextREDcond = fextRED_d -Hdr*(Hrr\fextRED_r) ; 
% 4) Dirichlet condition vector (condensed)
  %   { {\cCOND{} \defeq \cBC{}_d -  \Hqr{T}_{dr} \Hqr{-T}_{rr}\cBC{}_r}}
cCOND = cBC_d -Hrd'*(Hrr'\cBC_r) ;
% eINT (Rigid body forces at interfaces)
%  { {\eCONDint{} \defeq \Tint{}_r \Hqr{-1}_{rr} \fextRED{}_r}}
eCONDint = Tr*(Hrr\fextRED_r) ; 
%%%%%%
% Term bR, bP  
% \bMAT{r} = {\Hcond{T}}\KdomRED{-1}  \fextREDcond{}{} - \cCOND{} 
bR = Hcond'*(KdomREDglo\fextREDcond) -cCOND ; 
bP = -eCONDint ; 
%% INTERFACE FORCES 
% ------------------
%  \fINTF \defeq  \Tcond{}\Amat{rr}^{-1} \bMAT{r} -  \bMAT{p}  
fINTF = Tcond*(Arr\bR)-bP ; 

%%%% A) Solving for interface displacements 
pINT = Kintf\fINTF ;
%%%   B) Amplitude of seslf-equilibrated modes
% \rDEF{} = \Amat{rr}^{-1} \Par{\bMAT{r}   -  \Tcond{T}\pINT{}}
rDEF = Arr\(bR-Tcond'*pINT) ; 
%%%% C) Deformational displ. amplitudes
% \qDEF{} = \KdomRED{-1} (\fextREDcond - \Hcond{} \rDEF{})
 qDEF = KdomREDglo\(fextREDcond-Hcond*rDEF) ; 
%%% D) RIGID BODY MODEs
%\qRB{} = \Hqr{-T}_{rr} \Par{\cBC{}_r  -   {\Hqr{T}_{dr}}\qDEF{} +   {\Tint{T}_r} \pINT{}}

qRB = Hrr'\(cBC_r-Hdr'*qDEF +Tr'*pINT); 
% E) REaction resultants 
% { { \rRB{} =   {\Hqr{-1}_{rr}} \Par{\fextRED{}_r -   {\Hqr{}_{rd}}\rDEF{}}} }
rRB =  Hrr\(fextRED_r-Hrd*rDEF); 


% 
% 
% 
% nROW(1) = size(KdomREDglo,1) ;
% nROW(2) =  size(Hqr_r,1)  ;
% nROW(3) =   size(Hqr_r,2)   ;
% if isempty(Tint)
%     nROW(4)  = 0 ;
% else
%     nROW(4) =   size(Tint,2)   ;
%     
% end
% 
% 
% nCOL(1) = size(KdomREDglo,2) ;
% nCOL(2) =  size(Hqr_r,1)  ;
% nCOL(3) =   size(Hqr_d,2)   ;
% if isempty(Tint)
%     nCOL(4) = 0   ;
% else
%     nCOL(4) =   size(Tint,2)   ;
%     
% end
% 
% nrows = sum(nROW) ;
% ncols = sum(nCOL) ;
% 
% A  =sparse(nrows,ncols) ;
% b = zeros(nrows,1) ;
% 
% DATAIN = DefaultField(DATAIN,'STABILIZATION_PARAMETER',0) ;
% 
% maxK = max(max(KdomREDglo)); 
% 
% gSTB = (DATAIN.STABILIZATION_PARAMETER)*maxK ;
% PcompRR = gSTB*Pcomp(INDrig,INDrig) ;
% PcompRD = gSTB*Pcomp(INDrig,INDdef) ;
% PcompDD = gSTB*Pcomp(INDdef,INDdef) ;
% bCOMPr = gSTB*bCOMP(INDrig);
% bCOMPd = gSTB*bCOMP(INDdef) ;
% Tall =  gSTB*Tall ; 
% 
% %%%% BLOCK 1,1
% irow =1 ; icol = 1;
% indROW =  1:sum(nROW(1:irow)) ;
% indCOL = 1:sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = KdomREDglo + PcompDD  ;
% b(indROW) = fextRED_d + bCOMPd ;
% %%%% BLOCK 1,3
% irow =1 ; icol = 3;
% indROW = 1:sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = Hqr_d ;
% %%%% BLOCK 1,2
% irow =1 ; icol = 2;
% indROW = 1:sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = PcompRD' ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%% BLOCK 2,1
% irow =2 ; icol = 1;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = PcompRD ;
% b(indROW) = bCOMPr ;
% 
%  %%%% BLOCK 2,2
% irow =2 ; icol = 2;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = PcompRR ;
% %%%% BLOCK 2,3
% irow =2 ; icol = 3;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = Hqr_r ;
% b(indROW) = fextRED_r ;
% 
% %%%% BLOCK 3,1
% irow =3 ; icol = 1;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = Hqr_d' ;
% %%%% BLOCK 3,2
% irow =3 ; icol = 2;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = Hqr_r' ;
% %%%% BLOCK 3,3
% irow =3 ; icol = 3;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
% A(indROW,indCOL) = Tall ;
% %%%% BLOCK 3,4
% if ~isempty(Tint)
%     irow =3 ; icol = 4;
%     indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
%     indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
%     A(indROW,indCOL) = Tint ;
% end
% 
% 
% %%%% BLOCK 4,3
% if ~isempty(Tint)
%     irow =4 ; icol = 3;
%     indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
%     indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
%     A(indROW,indCOL) = Tint' ;
% end
% 
% % BLOCK q,qbar,r 
%  irow =3 ; icol = 3;
%   indROW = 1:sum(nROW(1:irow)) ;
%     indCOL = 1:sum(nCOL(1:icol)) ;
%     Add =    A(indROW,indCOL)   ;
%     % Block T 
%  irow =3 ; icol = 4;
%   indROW =1:sum(nROW(1:irow)) ;
%     indCOL = (sum(nCOL(1:icol-1))+1):sum(nCOL(1:icol)) ;
%     Adp =    A(indROW,indCOL)   ;    
% 
% 
% x = A\b ;
% 
% irow =1 ;
% qDEF = x(1:nROW(irow)) ;
% irow =2 ;
% indROW = (sum(nCOL(1:irow-1))+1):sum(nROW(1:irow)) ;
% qRB = x(1:nROW(irow)) ;
% irow =3 ;
% rDOM = x(1:nROW(irow)) ;
% 
% 
