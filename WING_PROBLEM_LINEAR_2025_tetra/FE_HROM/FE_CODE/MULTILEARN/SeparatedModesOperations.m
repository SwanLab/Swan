function [COLUMNS_RVEloc,nTYPErve,COLUMNS_RVE] = SeparatedModesOperations(DATAINM,nDOM) ; 

DATAINM = DefaultField(DATAINM,'SEPARATED_MODES',[]) ; % = TYPE_RVE ;
iproject = 1; 
DATAINM.SEPARATED_MODES = DefaultField(DATAINM.SEPARATED_MODES,'COLUMNS_RVETYPE_TEST',[]) ; % = TYPE_RVE ;
if ~isempty(DATAINM.SEPARATED_MODES.COLUMNS_RVETYPE_TEST) && DATAINM.SEPARATED_MODES.ACTIVE ==1
    nTYPErve = length(DATAINM.SEPARATED_MODES.COLUMNS_RVETYPE_TEST) ;
    COLUMNS_RVEloc = {ones(1,nDOM)} ;
    COLUMNS_RVE = cell(1,nTYPErve) ;
    iacumPROJ = 0 ;
    [COLUMNS_RVE,COLUMNS_RVEloc ]=...
        ColumnsSeparatedRVEs(DATAINM,nTYPErve,nDOM,COLUMNS_RVE,iproject,iacumPROJ,COLUMNS_RVEloc);
    COLUMNS_RVEloc = COLUMNS_RVEloc{1} ;
else
    COLUMNS_RVEloc = ones(1,nDOM) ;
    nTYPErve = 1 ;
    COLUMNS_RVE = {1:nDOM} ;
end