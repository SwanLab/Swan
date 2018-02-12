function L = update_lambda(problembsc,lambda,penalty,eta,constr,flag_up_date_lambda)
tipo = problembsc.lagrange_update;
tipo_augmented = problembsc.augmented_type;
L = lambda;

global lambda_max ;
global lambda_min ;
global contador_lamb;
global constr_old;

if flag_up_date_lambda
    switch tipo
        case 'AUGMENTED'
            %L = max(0,L + element.material.penalty*constr);
            
            switch tipo_augmented
                
                case 'NOVOTNY'
                    
                    
                    L = L + penalty*constr;
                    
                    
                    
                
                case 'BIS'
                    

                    if contador_lamb <= 1
                    
                        if constr  >= 0 
                            L = lambda*1.1;
                        else
                            L = lambda/1.1;
                        end
                        
                        if constr <= 0 && constr_old >= 0
                            if contador_lamb == 0
                                lambda_max = lambda;
                            end
                            contador_lamb = contador_lamb + 1;
                        elseif constr >= 0 && constr_old <= 0
                            if contador_lamb == 1
                                lambda_min = lambda;
                            end
                            contador_lamb = contador_lamb + 1;
                        end
                        
                    else
                        if contador_lamb >= 2
                            if constr <= 0
                                lambda_max = L;
                            else
                                lambda_min = L;
                            end
                        end
                        L = (lambda_max + lambda_min)/2;
                        contador_lamb = contador_lamb + 1;
%                         if constr <= 0
%                             lambda_max = L;
%                         else
%                             lambda_min = L;
%                         end
                    end

                       
                case 'THETA_DEPENDANCY'
                        L = L + penalty*constr;
                case 'NOCEDAL'
                    
                    if abs(constr) <= eta
                        L = L + penalty*constr;
                    else
                        L = L;
                    end
                    
                case 'INCREASE'
                    %L = L + 1/penalty*constr;
                    L = L + penalty*constr;
                    
            end
        case 'POTENCIAL'
    end
    
    
end
