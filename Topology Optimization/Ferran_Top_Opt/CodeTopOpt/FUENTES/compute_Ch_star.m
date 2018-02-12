function Ch_star = compute_Ch_star(matprop)
E_plus = matprop.E_plus;
E_minus = matprop.E_minus;
nu_plus = matprop.nu_plus;
nu_minus = matprop.nu_minus;

C_Cstar_case = 'nu_0_6'; % 'Seba' 'nu_0_6' 'negative_poisson' 'nu_0_8'

switch C_Cstar_case
    case 'negative_poisson'
        kappa_f = @(E,nu) E/2*(1-nu);
        mu_f = @(E,nu) E/2*(1-nu);
        
        k_plus = kappa_f(E_plus,nu_plus);
        mu_plus = mu_f(E_plus,nu_plus);
        
        k_minus = kappa_f(E_minus,nu_minus);
        mu_minus = mu_f(E_minus,nu_minus);

        
        nu = @(k,mu) (k-mu)/(k+mu);
        E = @(k,mu) (4*k*mu)/(k+mu);
        C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
        
        kappa_nu_min = k_minus;
        mu_nu_min = mu_plus;
        
        nu_min = nu(kappa_nu_min,mu_nu_min);
        E_nu_min = E(kappa_nu_min,mu_nu_min);
        C_nu_min = C(E_nu_min,nu_min);
        Ch_star = C_nu_min;
        
    case 'nu_0_6' %From Sigmund Thesis% rho = 0.38
        C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
        nu = -0.6;
        E = (1-nu*nu)*0.04;
        Ch_star = C(E,nu);
        
    case 'Seba' % Es=0.08; nus=-0.25
        
     Ch_star = [0.0853    -0.0213       0;
                -0.0213    0.0853       0;
                    0         0    0.0533];
                
    case 'nu_0_8' %From Sigmund Thesis% rho = 0.25
        C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
        nu = -0.8;
        E = (1-nu*nu)*0.02;
        Ch_star = C(E,nu);
       
end




end