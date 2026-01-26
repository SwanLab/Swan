function var_elem = AssignOneEIFelement_PROP(ngaus_element,PROPMATLOC,MaterialTypeLocal,NAMEVAR)


var_elem = zeros(ngaus_element,1) ;

for jmat = 1:length(PROPMATLOC)
    % Loop over number of materials of an EIF element
    var =PROPMATLOC(jmat).(NAMEVAR) ; %
    POINTSmat = find(MaterialTypeLocal == jmat) ;
    
    
    % INDEXES ROWS
 %   IND_POINTS = small2large(POINTSmat,nstrain) ;
    var_elem(POINTSmat,:)  = var ;
    
    
end