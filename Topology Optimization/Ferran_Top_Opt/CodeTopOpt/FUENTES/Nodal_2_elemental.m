function [gamma_gp,A_nod2gauss] = Nodal_2_elemental(gamma,element,dim,problembsc)
    if element.smoothing
        [gamma_gp,A_nod2gauss] = interpol(gamma,element,dim,problembsc);
    else
        gamma_gp = gamma;
        A_nod2gauss = eye(size(gamma));
    end
end
