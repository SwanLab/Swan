function CENTROID = ComputeCentroid(DATA3D,NODESfaces)
if nargin == 0
    load('tmp1.mat')
end

nface= length(NODESfaces) ;
CONNECTb_faces = cell(nface,1) ;
CENTROID =  cell(nface,1) ;
    ndim = size(DATA3D.COOR,2) ; 

for iface = 1:nface
    [dummy, setBelemLOC]= ElemBnd(DATA3D.CNb,NODESfaces{iface}); % elements face "iface"
    % Connectivities faces f1 and f2
    CONNECTb_faces{iface} = DATA3D.CNb(setBelemLOC,:) ;
    % The above is the connectivity matrix for the nodes of face "iface"
    COOR_FACE = DATA3D.COOR(NODESfaces{iface},:) ;
    [Nst,wST] = GeometricMassMatrixSurface(CONNECTb_faces{iface},NODESfaces{iface},COOR_FACE,DATA3D.TypeElementB) ;
    wSTdiag = CompWeightDiag(wST,1)  ;
    Mst = (wSTdiag*Nst)'*Nst ;
    % Recomputing centroid
    CentroidFA = zeros(1,ndim) ;
    AREA = sum(wST) ;
    for idim = 1:ndim
        CentroidFA(idim) = wST'*(Nst*COOR_FACE(:,idim))/AREA ;
    end
    CENTROID{iface} = CentroidFA ;
end