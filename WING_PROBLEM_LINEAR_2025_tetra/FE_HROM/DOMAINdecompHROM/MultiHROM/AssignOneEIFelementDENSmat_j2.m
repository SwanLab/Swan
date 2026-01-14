function dens_elem = AssignOneEIFelementDENSmat_j2(DATA,PROPMATLOC,MaterialTypeLocal)


dens_elem = zeros(DATA.MESH.ngaus_RHS,1) ;

for jmat = 1:length(PROPMATLOC)
    % Loop over number of materials of an EIF element
    dens =PROPMATLOC(jmat).dens ; %
    POINTSmat = find(MaterialTypeLocal == jmat) ;
    
    
    % INDEXES ROWS
 %   IND_POINTS = small2large(POINTSmat,nstrain) ;
    dens_elem(POINTSmat,:)  = dens ;
    
    
end