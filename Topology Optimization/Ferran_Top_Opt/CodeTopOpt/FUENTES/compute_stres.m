function stres = compute_stres(dim,strain,Ce,element,ptype)

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;
stres = zeros(nstre,nelem);

switch ptype
    case '2D'
        switch element.material.subtype
            case 'PLANESTRAIN'
                stres = zeros(nstre+1,nelem);
                epoiss = ones(1,nelem)*element.material.poiss;
                
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:) = stres(istre,:) + squeeze(Ce(istre,jstre,:))'.*strain(jstre,:);
                    end
                end
                stres(4,:) = epoiss.*(stres(1,:)+stres(2,:));

            case 'PLANESTRES'
                stres = zeros(nstre,nelem);
                
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:) = stres(istre,:) + squeeze(Ce(istre,jstre,:))'.*strain(jstre,:);
                    end
                end
                
         end
    case '3D'
        stres = zeros(nstre,nelem);
        
        for istre=1:nstre
            for jstre=1:nstre
                stres(istre,:) = stres(istre,:) + squeeze(Ce(istre,jstre,:))'.*strain(jstre,:);
            end
        end
end

end
