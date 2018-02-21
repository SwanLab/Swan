function stres = compute_stress_fourier(strain,Ce,dim)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

stres = zeros(nstre,nelem);
switch ptype
    case '2D'
                 stres = zeros(nstre,nelem);
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:) = stres(istre,:) - squeeze(Ce(istre,jstre,:))'.*strain(jstre,:);
                    end
                end
 
    case '3D'
end


end 