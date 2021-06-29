classdef ChTargetFactory < handle
    
    methods (Access = public, Static)
        
        function Ch_star = create(cParams)
            
            switch cParams.type
                case 'negative_poisson'
                    
                    E1  = cParams.E_plus;
                    E0  = cParams.E_minus;
                    nu1 = cParams.nu_plus;
                    nu0 = cParams.nu_minus;                    
                    kappa_f = @(E,nu) E/2*(1-nu);
                    mu_f = @(E,nu) E/2*(1-nu);
                    
                    k_plus = kappa_f(E1,nu1);
                    mu_plus = mu_f(E1,nu1);
                    
                    k_minus = kappa_f(E0,nu0);
                    mu_minus = mu_f(E0,nu0);
                    
                    
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
                    
                case 'Vfrac07'
                    Ch_star =[
                        0.4256    0.2837         0
                        0.2837    0.4256         0
                        0         0    0.1419];
                    
                case 'Vfrac06'
                    Ch_star =[
                        0.2909    0.1940         0
                        0.1940    0.2909         0
                        0         0    0.0970];
                    
                case 'Vfrac05'
                    Ch_star =[
                        0.1892    0.1261         0
                        0.1261    0.1892         0
                        0         0    0.0631];
                    
                case 'Vfrac04'
                    Ch_star =[
                        0.1141    0.0761         0
                        0.0761    0.1141         0
                        0         0    0.0380];
                    
                case 'Vfrac04b'
                      Ch_star = [ 0.778999640844672   0.237468158897016  -0.000000038110874
                                  0.237468158897016   0.778999640844672  -0.000000038110874
                                  -0.000000038110874  -0.000000038110874   0.234110861512206];
                    
                case 'Vfrac03'
                    Ch_star =[
                        0.0611    0.0407         0
                        0.0407    0.0611         0
                        0         0    0.0204];
                    
                case 'HorizontalRectangleInclusion'
                    Ch_star =[0.498384287901812   0.000858512870689   0.000000000000000
                              0.000858512870689   0.002575540644965   0.000000000000000
                              0.000000000000000   0.000000000000000   0.000862242916306 ];                    

                      
                    
                case 'Composite'
                    Ch_star =[1 0.15 0;
                        0.15 0.5 0;
                        0 0 0.2];
                case 'HoneyComb'
                    Ch_star =0.094*[1 0.75 0
                        0.75 1 0
                        0 0 0.125];
                case 'InvertedHoneyComb'
                    Ch_star =0.08*[1 -0.5 0
                        -0.5 1 0
                        0 0 0.06];
                case 'AcousticZeroShearA'
                    Ch_star =[1 1 0;
                        1 1 0;
                        0 0 0];
                case 'NegativePoiss06'
                    Ch_star =0.04*[1 -0.6 0;
                        -0.6 1 0;
                        0 0 0.8];
            end
        end
    end
    
end
