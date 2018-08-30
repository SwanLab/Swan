function phi0 = init_phi(dim,element,coordinates,Msmooth,problembsc)

phi_ini = -2*ones(dim.nelem,1);
if element.smoothing
    [coordinatesn,coordinatesa] = init_coord(coordinates);
    phi0 = smooth(dim.nelem,dim.npnod,dim.ndime,dim.nnode,phi_ini',Msmooth,element,coordinatesn,coordinatesa,problembsc);
    
else
    phi0(:,1) = phi_ini;
    
end
phi0(element.initial_holes) = 1;

end