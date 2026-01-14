function P = BoundaryIntefaceForcesLocal_kinconst(nlines,NODES_LINES,ndim,nnode,DATAROM,MESH2D,DOFsKEEP,FORCES)

if nargin == 0
    load('tmp2.mat')
end

% This function should be VECTORIZED
% -----------------------------------

P = zeros(nnode*ndim,1) ;
M = DATAROM{1}.MassMinterfaces ;  % Mass matrices each interface
BasisINTF = DATAROM{1}.BasisInt.BasisINTF  ; % V
BasisINTFall = DATAROM{1}.BasisIntRB  ; % Rigid body modes
DOFSface =  DATAROM{1}.BasisInt.DOFsFACE_V ;
[ndofNODES_U, ~]= cellfun(@size,BasisINTFall) ;
nRB= size(BasisINTFall{1},2) ;
nnodeE = length(BasisINTFall) ; 
%
iini = 1;
DOFS_NODE = cell(size(ndofNODES_U)) ;
for i= 1:length(ndofNODES_U)
    ifin = iini + ndofNODES_U(i)-1  ;
    DOFS_NODE{i} = iini:ifin ;
    iini = ifin +1 ; % ndofNODES_U{i} -1 ;
end

for iline = 1:nlines   % Loop over number of lines/faces defined in GID
    NODES = NODES_LINES{iline} ;  % (Midside) - Nodes defining the line/face
    ftracINTFlocRES =  FORCES.LINE{iline}(1:nRB) ; % Components of applied force (3 for 2D, 6 for 3D)
    
    for inodeMIDSIDE = 1:length(NODES)  % Loop over nodes defining the LINE/FACE
        nodeMID = NODES(inodeMIDSIDE) ; % Midside node under study.
        [ielem,iface] =  find(nodeMID == MESH2D.CN) ;
        % ---> \BasisINTF{e}{i} \Mintf{e}{i} \BasisINTFrbALL{e}{i}
        BasisINTF_i = BasisINTF(DOFS_NODE{iface},:) ;
        Mintf_i = M{iface} ;
        ftracCDOM_i = (BasisINTF_i'*Mintf_i*BasisINTFall{iface})*ftracINTFlocRES ;
        for anod=1:nnodeE
                a = DOFSface{nnodeE} ;
                Anod = MESH2D.CN(ielem,anod) ;  A = Nod2DOF(Anod,ndim) ;
                A = A(1:length(DOFSface{nnodeE})) ; 
                %%%%%
                P(A) = P(A) + ftracCDOM_i(a) ;
        end
        
        
        
    end
    
    
    
    
    
    
    %     DOFS = small2large(NODES_LINES{iline},ndim) ;
    %     FORCESloc = zeros(ndim,1) ;
    
    %
    %
    %     for e = 1:nelem
    %         elemtype = MESH2D.MaterialType(e) ;
    %         Ke = DATAROM{elemtype}.Kskel ;
    %         for anod=1:nnodeE
    %             a = Nod2DOF(anod,ndim) ;
    %             for bnod= 1:nnodeE
    %                 b = Nod2DOF(bnod,ndim) ;
    %                 Anod = MESH2D.CN(e,anod) ;  A = Nod2DOF(Anod,ndim) ;
    %                 Bnod = MESH2D.CN(e,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
    %                 %%%%%
    %                 K(A,B) = K(A,B) + Ke(a,b) ;
    %             end
    %         end
    %     end
    %
    %
    %     FORCESloc(1:length(FORCES.LINE{iline}))  =     FORCES.LINE{iline} ;
    %     nnodes = length(NODES_LINES{iline}) ;
    %     FORCES_rep = repmat(FORCESloc,nnodes,1) ;
    %     P(DOFS) = FORCES_rep ;
end


if ~isempty(DOFsKEEP)
    P = P(DOFsKEEP) ;
end