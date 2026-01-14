function celasglo_elem = AssignOneEIFelementELASTICITYmatCOROT(nstrain,DATA,PROPMATLOC,MaterialTypeLocal,typePROBLEM,ncecmINT_forces)


celasglo_elem = zeros(nstrain*length(MaterialTypeLocal),nstrain) ;

for jmat = 1:length(PROPMATLOC)
    % Loop over number of materials of an EIF element
    celas3D =PROPMATLOC(jmat).ElasticityMatrix ; %
    INVcelas3D = inv(celas3D) ;
    POINTSmat = find(MaterialTypeLocal == jmat) ;
    switch typePROBLEM
        case 'pstrain'
            if nstrain == 3
                rowcol = [1 2 6] ;
            elseif nstrain == 4
                rowcol = [1 2 6 3] ;
            else
                error('Option not implemented')
            end
            celas = celas3D(rowcol,rowcol) ;
            
        case 'pstress'
            if nstrain == 3
                rowcol = [1 2 6] ;
                celasINV3D = inv(celas3D) ;
                celasINV = celasINV3D(rowcol,rowcol) ;
                celas = inv(celasINV) ;
            else
                error('Option not implemented')
            end
        case '3D'
            celas = celas3D ;
    end
    
    % INDEXES ROWS
    IND_POINTS = small2large(POINTSmat,nstrain) ;
    celasglo_elem(IND_POINTS,:)  = repmat(celas,length(POINTSmat),1)  ;
    
    
end