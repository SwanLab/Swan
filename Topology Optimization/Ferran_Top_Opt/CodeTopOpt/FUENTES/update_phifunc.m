function [phifunct] = update_phifunc(theta,kappa,phifunct,g_nodal,norm_g,remesh,phiremesh,algorithm_update)

   switch remesh
        case 1
            phifunct = phiremesh;
        case 0
            
            switch algorithm_update
                case 'AMSTUTZ'
                    if theta < 1e-14
                        theta = 1e-12;
                    end
                    beta1 = sin((1-kappa)*theta)/sin(theta);
                    beta2 = sin(kappa*theta)/sin(theta);
                    
                    phifunct = beta1*phifunct + beta2*g_nodal/norm_g;
                    
                case 'BCN'
                    beta1 = 1;
                    beta2 = (cos(theta)+sqrt(cos(theta)^2-(kappa*cos(theta))^2));
                    
                    phifunct = 1/(1-kappa*cos(theta))*(beta1*phifunct + beta2*g_nodal/norm_g);
                    
                case 'BCN_gphi_SI'
                    beta1 = 1/(1+2*kappa*norm_g*cos(theta)+kappa^2*norm_g^2)^0.5;
                    beta2 = kappa/(1+2*kappa*norm_g*cos(theta)+kappa^2*norm_g^2)^0.5;
                    phifunct = beta1*phifunct + beta2*g_nodal;
                    
                    
                case 'BCN_gphi_EX'
                    txi = kappa*norm_g;
                    beta1 = -txi*cos(theta)+sqrt(1-txi^2*sin(theta)^2);
                    beta2 = kappa;
                    phifunct = beta1*phifunct + beta2*g_nodal;
                    
                    
                case 'BCN_phig_SI'
                    txi = 1/norm_g*(cos(theta)-sqrt((kappa-1-sin(theta))*(kappa-1+sin(theta))));
                    beta1 = 1/(1-kappa);
                    beta2 = -txi/(1-kappa);
                    phifunct = beta1*phifunct + beta2*g_nodal;
                    
                    
                case 'BCN_phig_EX'
                    beta1 = 1 + kappa;
                    beta2 = -((1+kappa)*cos(theta)-sqrt(1-((1+kappa)*sin(theta))^2))/norm_g;
                    phifunct = beta1*phifunct + beta2*g_nodal;
            end
     
            
            
            
    
    end
  
end

