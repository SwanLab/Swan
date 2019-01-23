function [theta,norm_dif_rel] = call_theta_level_set(phi,g,scalar_product,algorithm)

switch algorithm
    case 'level_set'
        norm_phi = sqrt(scalar_product(phi,phi));
        norm_g = sqrt(scalar_product(g,g));
        norm_g_f = sqrt(scalar_product(g/norm_g -phi,g/norm_g -phi));
        scl_phi_g = scalar_product(phi,g);
        theta = real(acos(scl_phi_g/(norm_phi*norm_g)));
        norm_dif_rel = norm_g_f;
    case 'Projected_gradient'
        theta = 0;
        norm_dif_rel = 1;
end

end