
function Fbe_glo = AddCellsValues(Fbe_glo,ELEMS_LOC,FA,ssAA)

FA = cell2mat(Fbe_glo(ELEMS_LOC)) + repmat(cell2mat(FA),length(ELEMS_LOC),1) ;
FA = mat2cell(FA,ssAA(1)*ones(length(ELEMS_LOC),1),1) ;
Fbe_glo(ELEMS_LOC) = FA ;

end
