function [stres,strain] = hooke_law(dim,Bmat,element,coordinatesn,coordinatesa,ptype,Ce)
% The most classical constitutive model, the hooke's law
% sigma = C:strain

strain = compute_strain(dim,Bmat,element,coordinatesn,coordinatesa);
stres = compute_stres(dim,strain,Ce,element,ptype);


end


