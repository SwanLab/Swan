function  [DATAwmethod]  = DeterminePcompBC(nrows,nDOM,COLloc,BasisUrb,BasisUdef,uBAR,alphaBC,f,nRB,nMODES_all)

  % Dimensions of matrix PcompBC
    PcompBC = sparse(nrows,nrows) ;
    bCOMPbc = zeros(nrows,1) ;
    iROWS = 0 ;
    mCOMPbc = 0 ; 
    for idom = 1:nDOM
        % Loop over domains (i)
        itype_i = COLloc(idom) ;  % Type of RVE
        BasisUdom_i = [BasisUrb BasisUdef{itype_i}] ;  % Basis matrix of displacements for this type of RVE
        nMODES_i = size(BasisUdom_i,2) ;  % Number of modes
        indROWS = iROWS + (1:nMODES_i) ; iROWS = iROWS + nMODES_i ; % Rows indices
        iCOLS = 0 ;
        %%
        % -----------
        indFACES = find(alphaBC(idom,:)==1) ;
        for iface = 1:length(indFACES)
            fLOC = cell2mat(f(indFACES(iface))) ;
            bCOMPbc(indROWS) = bCOMPbc(indROWS) + BasisUdom_i(fLOC,:)'*uBAR{idom,indFACES(iface)} ;
            mCOMPbc = mCOMPbc +uBAR{idom,indFACES(iface)}'*uBAR{idom,indFACES(iface)}  ;
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
                indFACES = find(alphaBC(idom,:)==1) ;
                fLOC = cell2mat(f(indFACES)) ;
                PcompBC(indROWS,indCOLS) = PcompBC(indROWS,indCOLS) + BasisUdom_i(fLOC,:)'*BasisUdom_i(fLOC,:) ;
            else
                
            end
        end
        
        %  INDrig =[INDrig;[iacum+(1:nRB)]'] ;
        iacum =  sum(nMODES_all(1:idom)) ;
    end
    
    DATAwmethod.PcompBC = PcompBC ; 
    DATAwmethod.bCOMPbc = bCOMPbc ; 
        DATAwmethod.mCOMPbc = mCOMPbc ; 
