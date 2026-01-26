function [COLUMNS_RVE,COLUMNS_RVEloc ]= ColumnsSeparatedRVEs(DATAIN,nTYPErve,size_dRVEloc,COLUMNS_RVE,iproject,iacumPROJ,COLUMNS_RVEloc)

if DATAIN.SEPARATED_MODES.ACTIVE == 1
    nSNAPS_loc = DATAIN.SEPARATED_MODES.COLUMNS{iproject} ;
    for itype = 1:nTYPErve
        INI =  nSNAPS_loc{itype}{1} ;
        FIN = nSNAPS_loc{itype}{2} ;
        if ischar(FIN) && strcmp(FIN,'end')
            FIN = size_dRVEloc ;
        elseif ischar(FIN) && strcmp(FIN,'end-1')
            FIN = size_dRVEloc-1 ;
        elseif isnumeric(FIN)
        else
            error('Option not available')
        end
        if ischar(INI) && strcmp(INI,'end')
            INI = size_dRVEloc ;
        elseif ischar(INI) && strcmp(INI,'end-1')
            INI = size_dRVEloc-1 ;
        elseif isnumeric(INI)
        else
            error('Option not available')
        end
        COLUMNS_RVE{itype} = [COLUMNS_RVE{itype},(INI:FIN)+iacumPROJ] ;
        COLUMNS_RVEloc{iproject}(INI:FIN) = itype;
    end
else
    COLUMNS_RVE{1} = [COLUMNS_RVE{1},length(COLUMNS_RVE{1})+(1:size_dRVEloc)] ; ; 
end