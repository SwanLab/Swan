classdef ShFunc_Chomog_fraction < ShFunc_Chomog
    properties (Access = private)
        alpha
        beta
    end
    methods
        function obj=ShFunc_Chomog_fraction(settings)
            obj@ShFunc_Chomog(settings);
            obj.alpha=settings.micro.alpha/norm(settings.micro.alpha);
            obj.beta=settings.micro.beta/norm(settings.micro.beta);
        end
        function computef(obj,x,physicalProblem,interpolation,filter)
            rho=filter.getP0fromP1(x);
            matProps=interpolation.computeMatProp(rho);
            physicalProblem.setMatProps(matProps);
            physicalProblem.computeChomog;
            obj.setPhysicalData(physicalProblem.variables);
            %Cost
            inv_matCh = inv(obj.Chomog);
            alpha_beta = obj.projection_Chomog(inv_matCh,obj.alpha,obj.beta);
            beta_alpha = obj.projection_Chomog(inv_matCh,obj.beta,obj.alpha);
            alpha_alpha = obj.projection_Chomog(inv_matCh,obj.alpha,obj.alpha);
            beta_beta = obj.projection_Chomog(inv_matCh,obj.beta,obj.beta);
            costfunc = alpha_beta/alpha_alpha + beta_alpha/beta_beta;
            
            %Gradient
            obj.compute_Chomog_Derivatives(physicalProblem.dim.nstre,physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus,x,interpolation,filter);
            
            beta1 = alpha_alpha*obj.beta - alpha_beta*obj.alpha;
            beta2 = beta_beta*obj.alpha - beta_alpha*obj.beta;
            g1 = obj.derivative_projection_Chomog(inv_matCh,obj.alpha,beta1,obj.Chomog_Derivatives,physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus,physicalProblem.dim.nstre);
            g2 = obj.derivative_projection_Chomog(inv_matCh,obj.beta,beta2,obj.Chomog_Derivatives,physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus,physicalProblem.dim.nstre);
            gradient = g1/(alpha_alpha)^2 + g2/(beta_beta)^2;
            
            mass=filter.Msmooth;
            gradient=filter.getP1fromP0(gradient(:));
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