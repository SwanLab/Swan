function  [P  bCOMP mCOMP INDrig INDdef DATAwmethod] = ...
    AssemblyPcompTILING2D3D(BasisUrb,BasisUdef,NODESfaces,uBAR,alphaBC,COL,COLloc,ndim,...
    THETAfaces,GAMMAfacesALL,DATAINM) ;
% See Implementation.pdf
%dbstop('6')
if nargin == 0
    load('tmp2.mat')
end
GAMMAfaces = GAMMAfacesALL.FACES ;
nRB = size(BasisUrb,2) ; % Number of rigid body modes
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nDEFrve ]= cellfun(@size,BasisUdef) ;  % Number of deformation modes
nMODESrve = nDEFrve + nRB ; % Total number of modes
nrows  = sum(nDOMrve.*nMODESrve) ;
P = sparse(nrows,nrows) ;
bCOMP = zeros(nrows,1) ;
mCOMP = 0 ;
INDrig = [] ;
iacum = 0 ;
nDOM = size(uBAR,1) ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisUdef)
    nMODES_all(COL{itype}) = size(BasisUdef{itype},2) + nRB ;
end
f =cell(size(NODESfaces));
for iface =1:length(NODESfaces)
    f{iface} = small2large(NODESfaces{iface},ndim) ;
end

ALPHA_THETA = alphaBC + THETAfaces ;
INDrig =[] ;
iROWS =0 ;

%%% Energy minimizatio
% if ~isempty(Hbar)
%             BasisUrb = Hbar*BasisUrb ;
% 
%     for iii  =1:length(BasisUdef)
%         BasisUdef{iii} = Hbar*BasisUdef{iii} ;
%     end
%     
%     for idom = 1:size(uBAR,1)
%         for jface = 1:size(uBAR,2)
%             if ~isempty(uBAR{idom,jface})
%                 fLOC =  (f{jface}) ;
%                 uBAR{idom,jface} =  Hbar(fLOC,fLOC)*uBAR{idom,jface};
%             end
%         end
%     end
%     
% end


for idom = 1:nDOM
    % Loop over domains (i)
    itype_i = COLloc(idom) ;  % Type of RVE
    BasisUdom_i = [BasisUrb BasisUdef{itype_i}] ;  % Basis matrix of displacements for this type of RVE
    nMODES_i = size(BasisUdom_i,2) ;  % Number of modes
    indROWS = iROWS + (1:nMODES_i) ; iROWS = iROWS + nMODES_i ; % Rows indices
    iCOLS = 0 ;
    %%
    % Vector bCOMP and mCOMP
    % -----------
    indFACES = find(alphaBC(idom,:)==1) ;
    for iface = 1:length(indFACES)
        fLOC = cell2mat(f(indFACES(iface))) ;
        bCOMP(indROWS) = bCOMP(indROWS) + BasisUdom_i(fLOC,:)'*uBAR{idom,indFACES(iface)} ;
        mCOMP = mCOMP +uBAR{idom,indFACES(iface)}'*uBAR{idom,indFACES(iface)}  ;
    end
    %%% Matrix Pcomp
    for jdom = 1:nDOM
        % Loop over domains (j)
        itype_j = COLloc(jdom) ;  % Type of RVE
        nMODES_j = size(BasisUdef{itype_j},2) + nRB ; % Number of modes
        indCOLS = iCOLS + (1:nMODES_j) ; iCOLS = iCOLS + nMODES_j ;  % Columns indices
        if idom == jdom
            % Diagonal terms
            %    dbstop('58')
            indFACES = find(ALPHA_THETA(idom,:)==1) ;
            fLOC = cell2mat(f(indFACES)) ;
            P(indROWS,indCOLS) = P(indROWS,indCOLS) + BasisUdom_i(fLOC,:)'*BasisUdom_i(fLOC,:) ;
        else
            % Off-diagonal terms
            if GAMMAfaces(idom,jdom) > 0
                BasisUdom_j = [BasisUrb BasisUdef{itype_j}] ; % Basis matrix of displacements for this type of RVE
                fi = f{GAMMAfaces(idom,jdom)} ;
                fj = f{GAMMAfaces(jdom,idom)} ;
                P(indROWS,indCOLS) = P(indROWS,indCOLS) - BasisUdom_i(fi,:)'*BasisUdom_i(fj,:) ;
            end
            
        end
    end
    
    INDrig =[INDrig;[iacum+(1:nRB)]'] ;
    iacum =  sum(nMODES_all(1:idom)) ;
end

%dbstop('75')
INDdef = (1:nrows)';
INDdef(INDrig) =  [] ;
% Prr = P(INDrig,INDrig) ;
% Prd = P(INDrig,INDdef) ;
% Pdd = P(INDdef,INDdef) ;
% bCOMPr = bCOMP(INDrig);
% bCOMPd = bCOMP(INDdef) ;




%%%% Approach based on minimization of work of the faces of the domains
DATAwmethod =  [];
if DATAINM.MinimizationBoundaryWork > 0
    % Checking whether it is a 1D -tiling problem
    FACES_CONTACT = sum(THETAfaces,1) ;
    IND_CONTACT= find(FACES_CONTACT) ;
    if sum(IND_CONTACT-[1 3]) ~=0
        error('Option not compatible. ONly 1D-tiling problems are allowed (along x-axis)')
    end
    
    [DATAwmethod]  = DeterminePcompBC(nrows,nDOM,COLloc,BasisUrb,BasisUdef,uBAR,alphaBC,f,nRB,nMODES_all) ;
    
    
    
end

end

