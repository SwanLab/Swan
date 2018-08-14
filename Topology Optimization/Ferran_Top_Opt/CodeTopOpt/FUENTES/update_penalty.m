function pen = update_penalty(problembsc,penalty,constr,flag_up_date_lambda,vol_n,element)
tipo = problembsc.lagrange_update;
tipo_augmented = problembsc.augmented_type;
factor = problembsc.penalty_factor;
penalty_max = problembsc.penalty_max;

if flag_up_date_lambda
    switch tipo
        case 'AUGMENTED'
            %L = max(0,L + element.material.penalty*constr);
            
            switch tipo_augmented
                
                case 'NOVOTNY' 
                    pen = penalty;
                    
                case 'THETA_DEPENDANCY'
                    
                    %%%%%%%%%%%%%
                    %alpha = 0.1;
                    %pen = ((1-alpha)*penalty + alpha*((1-theta_n/pi)/(theta_n/pi))^(1));
                    
                    %pen = ((1-theta_n/pi)/(theta_n/pi))^(1);
                    
                    %%%%%%%%%%%%%%%%%%%%
%                     alpha = 5;
%                     pen_min = 1e-12;%0.00390625;
%                     pen_max = 5;
                    
                    %theta_max = 180/180*pi;
                    %theta_min = 1/180*pi;
                    
                    %delta_theta_n = theta_max - theta_n;
                    %delta_theta = theta_max - theta_min;
                    
%                     exp_theta_n = exp(alpha*delta_theta_n);
%                     exp_theta = exp(alpha*delta_theta);
%                     pen = pen_min*((1-exp_theta_n)/(1-exp_theta)*(pen_max/pen_min - exp_theta) + exp_theta_n);
                    
                    %element.material.Vfrac = 0.6;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    alpha = 5;
                    pen_min = 1e-12;%0.00390625;
                    pen_max = 5;
                    vol_max = max(1-element.material.Vfrac,element.material.Vfrac);
                    vol_min = 1e-3*element.material.Vfrac;
                    
                    
                    delta_theta_n = (vol_max - abs(vol_n-element.material.Vfrac));
                    delta_theta = (vol_max - vol_min);
                    
                    exp_theta_n = exp(alpha*delta_theta_n);
                    exp_theta = exp(alpha*delta_theta);
                    pen = pen_min*((1-exp_theta_n)/(1-exp_theta)*(pen_max/pen_min - exp_theta) + exp_theta_n);

                    

                  %%%%%%%%%%%%%%%%%%%%%
%                   if length(theta_t) > 1
%                     alpha = 0;    
%                     delta_theta = theta_t(end) - theta_t(end-1);
%                     pen = penalty - alpha*delta_theta;
%                     else 
%                     pen = penalty;
%                     end

                    
                case 'INCREASE'
                    pen = min(penalty_max,penalty/factor);
                    
                case 'BIS'
                    pen = 0;
                    
            end
        case 'POTENCIAL'
            
            
    end
else
    pen = penalty;
    
end
