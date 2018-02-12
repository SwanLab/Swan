function     [stres,strain] = compute_stres_strain(dim,Bmat,element,problembsc,coordinatesn,coordinatesa,Ce,d_u)

ftype = problembsc.phisical_type;
ptype = problembsc.problemtype;

    switch ftype
        case {'ELASTIC'}
            [stres,strain] = hooke_law(dim,Bmat,element, coordinatesn,coordinatesa,ptype,Ce);
        case {'THERMAL'}
            [stres,strain] = fourier_law(dim,element,Bmat,ptype,Ce_plus,Ce);
    end
    
end