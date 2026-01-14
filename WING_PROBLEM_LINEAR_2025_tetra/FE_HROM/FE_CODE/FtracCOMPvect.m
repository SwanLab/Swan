function Ftrac = FtracCOMPvect(COOR,CNb,TypeElementB,Fnod,Tnod,CONNECTb) ;
% This subroutine   returns the  external forces due to boundary tractions
%(Ftrac). Inputs: COOR, CNbGLO: cell array with the connectivity matrices for
%the boundary elements along each spatial direction.
%Fnod (nnode*ndim x 1): Point loads vector. tracBglo:
%Cell array containing the  matrices   with the values of the prescribed tractions
%for each spatial direction.  TypeElementB: Type of boundary element (linear for 2D).
%dbstop('9')
if nargin==0
    load('tmp.mat')
end



ndim = size(COOR,2); nnode = size(COOR,1) ;
Fdis = zeros(ndim*nnode,1) ;

% Cell 2 mat for connectivity matrix 
CONNECTb = ConnectivityMatrixCell2Mat(CONNECTb) ; 



for idim = 1: ndim
    if ~isempty(CNb) && ~isempty(CNb{idim})
        Fdis_i = FdisCOMPvect(COOR,CNb{idim},TypeElementB,Tnod{idim},idim,CONNECTb);
        Fdis = Fdis + Fdis_i ;
    end
end
Ftrac = Fnod + Fdis ;
end