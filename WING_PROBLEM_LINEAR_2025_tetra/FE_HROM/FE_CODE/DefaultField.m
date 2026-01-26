function  DATAGEN = DefaultField(DATAGEN,FIELDVAR,dEFval) ; 
% DefaultField adds the field FIELDVAR to structure array DATAGEN and
% assigns the specified value dEFval. 
if ~isfield(DATAGEN,FIELDVAR)
    DATAGEN.(FIELDVAR) = dEFval ; 
end