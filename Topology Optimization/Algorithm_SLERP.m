classdef Algorithm_SLERP < handle
    properties
        stop_Criteria_ls=1;
        stop_Criteria_opt=1;
        lambda=0;
        kappa=1;
        kfrac=2;
        max_constr_change=+Inf;
        kappa_min=1.0000e-15;
        penalty = 1;
        constr_tol=1;
        optimality_tol=0.0175;
        Msmooth
        Stiff_smooth
        epsilon_scalar_product_P1=0.03;
        fhtri
    end 
    methods
        function x=updateX(obj,x_ini,cost,constraint, physProblem, interpolation,filter)  
            obj.Msmooth=physProblem.computeMass(2);
            obj.Stiff_smooth=physProblem.computeKsmooth;
            cost.h_C_0=cost.value;
            iter=0;
            

            while(obj.stop_Criteria_opt)
                iter=iter+1;
                obj.plotX(x_ini,physProblem)
                volume = constraint.value;
                cost.computef(x_ini,physProblem,interpolation,filter);
                constraint.computef(x_ini,physProblem,interpolation,filter);
                obj.lambda = obj.lambda+obj.penalty*constraint.value;
                
                cost_ini = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                gradient_ini = constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
                
                theta = obj.computeTheta(x_ini,gradient_ini);
                while(obj.stop_Criteria_ls) 
                    
                    x_ls=obj.designVariableUpdate(x_ini,obj.kappa,theta,gradient_ini);  
                    
                    obj.plotX(x_ls,physProblem)
                    cost.computef(x_ls,physProblem,interpolation,filter);
                    constraint.computef(x_ls,physProblem,interpolation,filter);
                    
                    cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                    gradient_ls = constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;                    
                    theta_ls =obj.computeTheta(x_ls,gradient_ls);
                    volume_ls = constraint.value;
                    
                    incr_vol_ls = abs(volume_ls - volume);
                    incr_cost = (cost_ls - cost_ini)/abs(cost_ini);
                    
                    obj.kappa = obj.kappa/obj.kfrac;
                    obj.stop_Criteria_ls = ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.kappa <= obj.kappa_min);              
                end
                obj.stop_Criteria_ls=1;
                x_ini=x_ls;
                obj.kappa=1;
                active_constr = obj.penalty > 0;
                obj.stop_Criteria_opt = theta >= obj.optimality_tol || any(abs(constraint.value(active_constr)) > obj.constr_tol(active_constr));
            end
            x=x_ini;
        
        end
        function theta=computeTheta(obj,phi,g)
            norm_phi = sqrt(obj.scalar_product(phi,phi));
            norm_g = sqrt(obj.scalar_product(g,g));
            %norm_g_f = sqrt(scalar_product(g/norm_g -phi,g/norm_g -phi));
            scl_phi_g = obj.scalar_product(phi,g);
            theta = real(acos(scl_phi_g/(norm_phi*norm_g)));
            %norm_dif_rel = norm_g_f;
        end
        function phi=designVariableUpdate(obj,design_variable,kappa,theta,gradient)
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product(gradient,gradient));
            
            beta1 = sin((1-kappa)*theta)/sin(theta);
            beta2 = sin(kappa*theta)/sin(theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
        end
        function sp=scalar_product(obj,f,g)
            sp=f'*(((obj.epsilon_scalar_product_P1)^2)*obj.Stiff_smooth+obj.Msmooth)*g;
        end
    end
    methods
        function plotX(obj,x,physicalProblem)
            
                gamma_nodal=x<0;
            if isempty(obj.fhtri)
                fh = figure;
                mp = get(0, 'MonitorPositions');
                select_screen=1;
                if size(mp,1) < select_screen
                    select_screen = size(mp,1);
                end
                width = mp(1,3);
                height = mp(1,4);
                size_screen_offset = round([0.7*width,0.52*height,-0.71*width,-0.611*height],0);
                set(fh,'Position',mp(select_screen,:) + size_screen_offset);
                obj.fhtri = trisurf(physicalProblem.mesh.connec,physicalProblem.mesh.coord(:,1),physicalProblem.mesh.coord(:,2),double(gamma_nodal), ...
                    'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                view([0,90]);
                colormap(flipud(gray));
                set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
            else
                set(obj.fhtri,'FaceVertexCData',double(gamma_nodal));
                set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
                drawnow;
            end  
        end
    end
    
end
