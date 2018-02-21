function [fext,fext_adjoint,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass,Group] = call_read_data(file_name,method,matprop,TYPE,perimeter_case)
[fext,fext_adjoint,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass,Group] = read_data_problem(file_name,TYPE,perimeter_case);

if exist('matprop','var') % set material properties
    element.E_plus = matprop.E_plus;
    element.E_minus = matprop.E_minus;
    element.nu_plus = matprop.nu_plus;
    element.nu_minus = matprop.nu_minus;
end


switch method
    case 'SIMP'
        pSIMPfactor_fun = @(nu_plus) max(2/(1-nu_plus),4/(1+nu_plus));
        pSIMPfactor = pSIMPfactor_fun(element.nu_plus);
        
        [element.mu_func,mu] = mu_interp_simp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus,pSIMPfactor);
        [element.kappa_func,kappa] = K_interp_simp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus,pSIMPfactor);
                
        [~,P] = polarization_interp(mu,kappa);
        element.polarization_sym_part_fourth_order_tensor = matlabFunction(P(3,3));
        element.polarization_doble_product_second_identity_tensor = matlabFunction(P(1,2));
        
        % simp-all for comparison
        [element.mu_func_simp_all,mu_simp_all] = mu_interp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus);
        [element.kappa_func_simp_all,kappa_simp_all] = kappa_interp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus);

        [~,P_simp_all] = polarization_interp(mu_simp_all,kappa_simp_all);
        element.polarization_sym_part_fourth_order_tensor_simp_all = matlabFunction(P_simp_all(3,3));
        element.polarization_doble_product_second_identity_tensor_simp_all = matlabFunction(P_simp_all(1,2));
        
    case 'SIMP_ALL'
        [element.mu_func,mu] = mu_interp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus);
        [element.kappa_func,kappa] = kappa_interp(element.gamma_minus,element.gamma_plus,element.E_minus,element.E_plus,element.nu_minus,element.nu_plus);
        
        [~,P] = polarization_interp(mu,kappa);
        element.polarization_sym_part_fourth_order_tensor = matlabFunction(P(3,3));
        element.polarization_doble_product_second_identity_tensor = matlabFunction(P(1,2));

    otherwise
        error('Method not valid.');
end

end