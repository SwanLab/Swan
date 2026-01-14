function  [Tdd, c, s,INDrig,INDdef,T]=  AssemblyTreacRESULTANT(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COL,COLloc,...
    THETAfaces,GAMMAentities,ndim,alphaBC,r,BasisRrb,COORref)

%dbstop('5')
if nargin == 0
    load('tmp5.mat')
    %GAMMAentities = GAMMAfaces ;
end

% if ~isempty(HbarINV)
%             BasisRrb = HbarINV*BasisRrb ;
% 
%     for iii  =1:length(BasisRdef)
%         BasisRdef{iii} = HbarINV'*BasisRdef{iii} ;
%     end
% end
BasisUrbLOC ={} ; 
for iface = 1:length(NODESbound.PLANE)
  %  warning('Not implemented yet !!!' )
    FACEnodes = NODESbound.PLANE{iface} ; % Nodes involved face
    COOR_FACE = COORref(FACEnodes,:) ;   
    COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity     
    COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
    BasisUrbLOC{iface} = ConstructBasisRigidBody(COORrel) ; %
    DOFSloc = small2large(FACEnodes,ndim) ; 
    
end

nRB= size(BasisRrb,2) ;
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nMODES ]= cellfun(@size,BasisRdef) ;  % Number of reaction modes
nMODES = nMODES + nRB ;
nrows  = sum(nDOMrve.*nMODES) ;
T = sparse(nrows,nrows) ;
%c = zeros(nrows,1);
%s = 0 ;
iROWS = 0 ;
iCOLS = 0 ;
nDOM = size(betaBC,1) ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) + nRB ;
end
INDrig = [] ;
iacum = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BETA_THETA = betaBC + THETAfaces ;
BT = cell(size(BETA_THETA)) ; % While BETA_THETA(idom,jface) = 1, if jface of idom is a Neuman boundary,
%   BT{idom,jface} contains the set of nodes
%   which are actually free of prescribed
%                               %   displacements
NODESfaces = NODESbound.PLANE ;
GAMMAfaces = GAMMAentities.FACES ;

for idom = 1:size(BT,1)
    for iface = 1:size(BT,2)
        if BETA_THETA(idom,iface) == 1
            nodesF = NODESfaces{iface} ; %
            % Checking wether some of nodes "nodesF" are Dirichlet nodes
            alphaLOC = alphaBC(idom,:) ;
            faceDIRIC = find(alphaLOC==1) ;
            for icheck = 1:length(faceDIRIC)
                ifaceREMOVE = faceDIRIC(icheck);
                nodesREMOVE = NODESfaces{ifaceREMOVE} ;
                [ NodesIntersect iA iB]= intersect(nodesF,nodesREMOVE) ;
                nodesF(iA) = [] ;
            end
            BT{idom,iface} = nodesF' ;
        end
    end
end




% -------------------------------------------

%
% f =cell(size(NODESfaces));
% for iface =1:length(NODESfaces)
%     f{iface} = small2large(NODESfaces{iface},ndim) ;
% end
% BTorig = BT;
% for idom =1:size(BT,1)
%     for iface = 1:size(BT,2)
%         BT{idom,iface} = small2large(BT{idom,iface},ndim)' ;
%     end
% end


INDrig =[] ;
iROWS =0 ;
for idom = 1:nDOM
    % Loop over domains (i)
    itype_i = COLloc(idom) ;  % Type of RVE
    BasisRdom_i =  [BasisRrb, BasisRdef{itype_i} ];  % Basis matrix of displacements for this type of RVE
    nMODES_i = size(BasisRdom_i,2) ;  % Number of modes
    indROWS = iROWS + (1:nMODES_i) ; iROWS = iROWS + nMODES_i ; % Rows indices
    iCOLS = 0 ;
    %%
    %     % Vector c , and scalar s% --------------------
    %     indFACES = find(betaBC(idom,:)==1) ;
    %     for iface = 1:length(indFACES)
    %         fLOC = cell2mat(f(indFACES(iface))) ;
    %         c(indROWS) = c(indROWS) + BasisRdom_i(fLOC,:)'*reactDOMrbGLO{idom}(fLOC) ;
    %         s = s + norm(reactDOMrbGLO{idom}(fLOC))^2 ;
    %     end
    %%% Matrix Pcomp, vector c
    for jdom = 1:nDOM
        % Loop over domains (j)
        itype_j = COLloc(jdom) ;  % Type of RVE
        nMODES_j = size(BasisRdef{itype_j},2) +nRB   ; % Number of modes
        indCOLS = iCOLS + (1:nMODES_j) ; iCOLS = iCOLS + nMODES_j ;  % Columns indices
        if idom == jdom
            % Diagonal terms
            indFACES = find(BETA_THETA(idom,:)==1) ;
            fLOC = cell2mat(BT(idom,indFACES)) ;
            fLOC = unique(fLOC) ; % Change to avoid repeated nodes
            fLOC = small2large(fLOC,ndim) ;
            T(indROWS,indCOLS) = T(indROWS,indCOLS) + BasisRdom_i(fLOC,:)'*BasisRdom_i(fLOC,:) ;
        else
            % Off-diagonal terms
            % FACES
            if GAMMAfaces(idom,jdom) > 0
                BasisRdom_j = [BasisRrb, BasisRdef{itype_j}] ; % Basis matrix of displacements for this type of RVE
                fi = BT{idom,GAMMAfaces(idom,jdom)} ;
                fj = BT{jdom,GAMMAfaces(jdom,idom)} ;
                fi = small2large(fi,ndim) ;
                fj = small2large(fj,ndim) ;
                T(indROWS,indCOLS) = T(indROWS,indCOLS) + BasisRdom_i(fi,:)'*BasisRdom_j(fj,:) ;
            end
            % Lines
            if GAMMAentities.LINES(idom,jdom) > 0
                BasisRdom_j = [BasisRrb, BasisRdef{itype_j}] ; % Basis matrix of displacements for this type of RVE
                indI = GAMMAentities.LINES(idom,jdom) ;
                indJ = GAMMAentities.LINES(jdom,idom) ;
                fi = NODESbound.LINES{indI} ;
                    fj = NODESbound.LINES{indJ} ;
                fi = small2large(fi,ndim) ;
                fj = small2large(fj,ndim) ;
                T(indROWS,indCOLS) = T(indROWS,indCOLS) + BasisRdom_i(fi,:)'*BasisRdom_j(fj,:) ;
            end
            % Points
            if GAMMAentities.POINTS(idom,jdom) > 0
                BasisRdom_j = [BasisRrb, BasisRdef{itype_j}] ; % Basis matrix of displacements for this type of RVE
                indI = GAMMAentities.POINTS(idom,jdom) ;
                indJ = GAMMAentities.POINTS(jdom,idom) ;
                fi = NODESbound.POINTS{indI} ;
                fj = NODESbound.POINTS{indJ} ;
                fi = small2large(fi,ndim) ;
                fj = small2large(fj,ndim) ;
                T(indROWS,indCOLS) = T(indROWS,indCOLS) + BasisRdom_i(fi,:)'*BasisRdom_j(fj,:) ;
            end
        end
    end
    INDrig =[INDrig;[iacum+(1:nRB)]'] ;
    iacum =  sum(nMODES_all(1:idom)) ;
    
end
INDdef = (1:nrows)';
INDdef(INDrig) =  [] ;
Trr = T(INDrig,INDrig) ;
Tdr = T(INDdef,INDrig) ;
Tdd = T(INDdef,INDdef) ;

 c = Tdr*r ;
s = r'*Trr*r ;


end

