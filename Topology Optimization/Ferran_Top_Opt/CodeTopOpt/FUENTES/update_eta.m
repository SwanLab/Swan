function eta = update_eta(problembsc,penalty,constr,eta,flag_up_date_lambda)
tipo = problembsc.lagrange_update;
tipo_augmented = problembsc.augmented_type;


if flag_up_date_lambda
    switch tipo
        case 'AUGMENTED'
            %L = max(0,L + element.material.penalty*constr);
            
            switch tipo_augmented
                case 'THETA_DEPENDANCY'
                    eta = eta;
                case 'NOCEDAL'
                    
                    if abs(constr) <= eta
                        eta = eta/((penalty)^(0.9));
                    else
                        eta = 1/((penalty)^(0.1));
                    end
                    
                case 'INCREASE'
                     eta = eta;
                    
            end
        case 'POTENCIAL'
    end
    
else
    eta = eta;
    
end
