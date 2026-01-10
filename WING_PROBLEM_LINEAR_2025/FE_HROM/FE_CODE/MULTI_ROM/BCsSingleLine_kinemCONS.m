function  [DOFr,dR] = BCsSingleLine_kinemCONS(NODES_LINES,ilines,DISP,ndim,BasisInt,MESH2D,nRB)

if nargin == 0 
    load('tmp2.mat')
    ilines = 3; 
end

%warning('This should be checkeD !!!!')

% Total number of DOFs  
nface = length(BasisInt.BasisINTFall_cell);

IndicesRB = BasisInt.IndicesRB ; % Indices DOFs of BasisINTFall
IndicesRB = mat2cell(IndicesRB,1,nRB*ones(1,nface));
DOFmFACE =  BasisInt.DOFsFACE   ;  % Master DOFs of BasisINTFall, face-wise
DOFmFACE_V =  BasisInt.DOFsFACE_V   ;  % Master DOFs of BasisINTF, face-wise
DISPLOC  =DISP.LINE{ilines} ; % Prescribed displacement ---rigid body type 3 or 6 entries. 

% First we have to determine which is the face associated to the line under
% study 
[~,iface]   = find(MESH2D.CN==NODES_LINES{ilines}(1)); 
IndicesRB = IndicesRB{iface} ;  % Indexes RB corresponding to this face

DOFmFACE =  BasisInt.DOFsFACE{iface} ;  % DOFS of BasisINTFall that are master DOFs
[DOFSrb IndicesRBloc]= intersect(DOFmFACE,IndicesRB) ;  % Find which of these DOFs are RIGID BODY dofs
nDOFsFACE = length(DOFmFACE)  ; % Number of total DOFs associated to the face under consideration 


DISPLOC_all  =cell(1,nDOFsFACE) ; 
DISPLOCrb  =DISP.LINE{ilines} ; % Data provided  by the user (Rigid body comp. 3 or 6)
DISPLOCrb = DISPLOCrb(IndicesRBloc) ; 

if any(~isempty(DISPLOCrb))
    % Some of the disp. prescribed by the user is not empty. Then we make
    % all DOFs zero 
    DISPLOC_all(:) = {0} ; 
    % And then set 
    DISPLOC_all(IndicesRBloc) = DISPLOCrb ; 
end

 DISPLOC = DISPLOC_all ; 

DOFS = small2large(NODES_LINES{ilines},ndim) ; % All DOFs involved line  (ndim --> nmax)
nnodes = length(NODES_LINES{ilines}) ;
% Known DOFs
dR = [] ;
DOFr = [] ;

ndimLOC = length(DISPLOC) ;
PRESCRIBED_ALL = 0 ;
for idim = 1:ndim
    if idim <=ndimLOC
        if ~isempty(DISPLOC{idim})
            DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
            dR = [dR ; DISPLOC{idim}*ones(nnodes,1)] ;
            PRESCRIBED_ALL = 1;
        end
    elseif idim > ndimLOC
        if PRESCRIBED_ALL == 1
            DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
            dR = [dR ; 0*ones(nnodes,1)] ;
        end
    end
end
DOFr = DOFS(DOFr) ;

end





%
%
% % RIGHT END
% %----------------
%
%
% NODE_2 = MESH2D.RIGHT_END_NODE ;
% DOFS_2= small2large(NODE_2,ndim) ;
% % Known DOFs
% r = [] ;
% DISPLOC  =DISP.RIGHT_END ;
% for idim = 1:length(DISPLOC)
%     if ~isempty(DISPLOC{idim})
%         r = [r; idim] ;
%         dR = [dR ; DISPLOC{idim}] ;
%     end
% end
% DOFr = [DOFr; DOFS_2(r)] ;
%
%
%
% DOFl = setdiff(1:ndim*nnode,DOFr) ;
%
%
