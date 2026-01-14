function IND_ELEM_TYPE = ChooseNumberDomainsToPrint(MESH3D,IND_ELEM_TYPE,iii)

if ~isempty(MESH3D.DOMAINSprint.INDICES)
        INDD =MESH3D.DOMAINSprint.INDICES{iii} ;
        switch  MESH3D.DOMAINSprint.TYPE_INDEX 
            case 'LOCAL'
                IND_ELEM_TYPE = IND_ELEM_TYPE(INDD) ; 
            case 'GLOBAL'
                IND_ELEM_TYPE = INDD ; 
            case 'NUMBER_OF_SLICES'
                INDDa = linspace(1,length(IND_ELEM_TYPE),INDD) ;
                INDDa = floor(INDDa) ;
                IND_ELEM_TYPE = IND_ELEM_TYPE(INDDa) ; 
        end    
    end