function Ftrac = FtracCOMP(COOR,CNb,TypeElementB,Fnod,Tnod) ;
% This subroutine   returns the  external forces due to boundary tractions
%(Ftrac). Inputs: COOR, CNbGLO: cell array with the connectivity matrices for
%the boundary elements along each spatial direction.
%Fnod (nnode*ndim x 1): Point loads vector. tracBglo:
%Cell array containing the  matrices   with the values of the prescribed tractions
%for each spatial direction.  TypeElementB: Type of boundary element (linear for 2D).
%dbstop('9')
if nargin==0
    load('tmp1.mat')
end
ndim = size(COOR,2); nnode = size(COOR,1) ;
Fdis = zeros(ndim*nnode,1) ;
for idim = 1: ndim
    if ~isempty(CNb{idim})
        Fdis_i = FdisCOMP(COOR,CNb{idim},TypeElementB,Tnod{idim},idim);
        Fdis = Fdis + Fdis_i ;
    end
end
Ftrac = Fnod + Fdis ;
end