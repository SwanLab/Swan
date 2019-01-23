function [stres_plus,stres_minus,strain] = fourier_law(dim,element,Bmat,ptype,Ce_plus,Ce_minus,d_u)
% The most classical constitutive model, the hooke's law
% sigma = C:strain

strain = compute_strain_fourier(dim,element,Bmat,d_u);
stres_plus = compute_stress_fourier(strain,Ce_plus,dim);
stres_minus = compute_stress_fourier(strain,Ce_minus,dim);


end



