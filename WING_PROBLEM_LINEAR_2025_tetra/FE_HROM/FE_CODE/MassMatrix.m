function M= MassMatrix(DATA,Nst,wSTs,ndim)
%dbstop('3')
if nargin == 0
    load('tmp.mat')
end
densGLO = repmat(DATA.densGLO',DATA.ngaus,1) ;
densGLO = densGLO(:) ;
densGLO = CompWeightDiag(densGLO,ndim) ;  ;
wDIAG = CompWeightDiag(wSTs,ndim)  ;
NstW = wDIAG*Nst ;
M = densGLO*Nst ;
M = NstW'*M ;