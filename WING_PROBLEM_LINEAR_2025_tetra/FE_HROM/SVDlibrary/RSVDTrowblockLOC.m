function [U,S,V] = RSVDTrowblockLOC(A,RELTOL,DATA)
% See comments in RSVDTrowblock.m 

if nargin == 0
    load('tmp.mat')
end
if nargin == 2
    DATA = [] ;
end

if length(RELTOL) == length(A)
    % Block-wise tolerance
    TOL_BLOCK = RELTOL ;
    TOL_GLO = 0 ;
else
    TOL_BLOCK = 0*ones(size(A)) ;
    TOL_GLO = RELTOL ;
end

DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ; % = 0;
[Q,L,gamma,NormA,ETIME] = RORTHrowblock(A,TOL_BLOCK,DATA) ;
% Q is an "exact" orthogonal basis matrix for the column space of A  y
% TOL_BLOCK =[0,0 ...0]; otherwise, it is an approximated basis matrix,
% with approximation given by TOL_BLOCK
NormA = norm(NormA) ;


e0 = TOL_GLO ;
eORTH = norm(gamma) ;
%if  size(Q,2) < ncolsTOT
%  if DATA2.RELATIVE_SVD ==1 && TOL_GLO >0 %
e0 = TOL_GLO*NormA;
% end
if e0 > eORTH;   e = sqrt(e0^2 - eORTH^2) ;   else   e = e0 ; end

%    dbstop('95')
DATA.RELATIVE_SVD = 0 ;
[sA,sB] = cellfun(@size,A,'UniformOutput',false) ;
nrowsTOT = sA{1} ;
ncolsTOT = sum(cell2mat(sB)) ;


DATA.SVD.MaxSize = max([nrowsTOT,ncolsTOT])  ;
[U,S,V,eSVDb] = SVDT(cell2mat(L),e,DATA);
U = Q*U ;
eSVD = sqrt(eORTH^2 + eSVDb^2) ;
%else
% A appears to be full rank
%   [U,S,V,eSVD] = SVDT(A,e0,DATA);

%end
