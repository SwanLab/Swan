function dens_elem = AssignOneEIFelementDENSmat(DATA,PROPMATLOC,MaterialTypeLocal)


dens_elem = zeros(DATA.MESH.ngaus_RHS,1) ;

for jmat = 1:length(PROPMATLOC)
    % Loop over number of materials of an EIF element
    dens =PROPMATLOC(jmat).Density ; %
    POINTSmat = find(MaterialTypeLocal == jmat) ;
    
    
    % INDEXES ROWS
 %   IND_POINTS = small2large(POINTSmat,nstrain) ;
    dens_elem(POINTSmat,:)  = dens ;
    
    
end


if any(dens_elem==0)
    error('Density cannot be equal to zero...Check function in which submaterials are defined')
end