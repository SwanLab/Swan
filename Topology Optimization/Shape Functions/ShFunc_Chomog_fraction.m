classdef ShFunc_Chomog_fraction < ShFunc_Chomog
    properties (Access = public)
        alpha
        beta
    end
    methods
        function obj=ShFunc_Chomog_fraction(settings)
            obj@ShFunc_Chomog(settings);
            obj.alpha=settings.micro.alpha/norm(settings.micro.alpha);
            obj.beta=settings.micro.beta/norm(settings.micro.beta);
        end
        function computef(obj,x)
            obj.rho=obj.filter.getP0fromP1(x);
            obj.matProps=obj.interpolation.computeMatProp(obj.rho);
            obj.physicalProblem.setMatProps(obj.matProps);
            obj.physicalProblem.computeChomog;
            obj.setPhysicalData(obj.physicalProblem.variables);
            
            %Cost
            inv_matCh = inv(obj.Chomog);
            alpha_beta = obj.projection_Chomog(inv_matCh,obj.alpha,obj.beta);
            beta_alpha = obj.projection_Chomog(inv_matCh,obj.beta,obj.alpha);
            alpha_alpha = obj.projection_Chomog(inv_matCh,obj.alpha,obj.alpha);
            beta_beta = obj.projection_Chomog(inv_matCh,obj.beta,obj.beta);
            costfunc = alpha_beta/alpha_alpha + beta_alpha/beta_beta;
            
            %Gradient
            obj.compute_Chomog_Derivatives(x);
            
            beta1 = alpha_alpha*obj.beta - alpha_beta*obj.alpha;
            beta2 = beta_beta*obj.alpha - beta_alpha*obj.beta;
            g1 = obj.derivative_projection_Chomog(inv_matCh,obj.alpha,beta1);
            g2 = obj.derivative_projection_Chomog(inv_matCh,obj.beta,beta2);
            gradient = g1/(alpha_alpha)^2 + g2/(beta_beta)^2;
            
            % !! NOT FROM FILTER !!
            mass=obj.filter.diffReacProb.element.M;
            
            gradient=obj.filter.getP1fromP0(gradient(:));
            gradient = mass*gradient;
            if isempty(obj.h_C_0)
                obj.h_C_0 = costfunc;
            end
            costfunc = costfunc/abs(obj.h_C_0);
            gradient=gradient/abs(obj.h_C_0);
%             obj.h_C_0 = costfunc;
            
            obj.value = costfunc;
            obj.gradient = gradient;
        end
    end
end