function [DOFr,DOFl,dR] = DirichletBNDCondBeam_mixed(DATAROM,MESH1D,DISP,DATAIN,ndim)

if nargin == 0
    load('tmp.mat')
end


nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

% LEFT END
% ----------
NODE_1 = MESH1D.LEFT_END_NODE ;
DOFS_1 = small2large(NODE_1,ndim) ;
% Known DOFs
dR = [] ;
r = [] ;
DISPLOC  =DISP.LEFT_END ;
for idim = 1:length(DISPLOC)
    if ~isempty(DISPLOC{idim})
        r = [r; idim] ;
        dR = [dR ; DISPLOC{idim}] ;
    end
end

DOFr = DOFS_1(r)  ;



% RIGHT END
%----------------


NODE_2 = MESH1D.RIGHT_END_NODE ;
DOFS_2= small2large(NODE_2,ndim) ;
% Known DOFs
r = [] ;
DISPLOC  =DISP.RIGHT_END ;
for idim = 1:length(DISPLOC)
    if ~isempty(DISPLOC{idim})
        r = [r; idim] ;
        dR = [dR ; DISPLOC{idim}] ;
    end
end
DOFr = [DOFr; DOFS_2(r)] ;



DOFl = setdiff(1:ndim*nnode,DOFr) ;

