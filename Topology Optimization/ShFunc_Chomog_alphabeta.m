classdef ShFunc_Chomog_alphabeta < ShFunc_Chomog
    properties (Access = private)
        alpha;
        beta;
    end
    methods
        function obj=ShFunc_Chomog_alphabeta(settings)
            obj@ShFunc_Chomog(settings);
            obj.alpha=settings.micro.alpha;
            obj.beta=settings.micro.beta;
        end
        function computef(obj,x,physicalProblem,interpolation,filter)
            %[obj.Chomog,obj.tstress,obj.tstrain] = physicalProblem.computeChomog;
            inv_matCh = inv(obj.Chomog);
            costfunc = obj.projection_Chomog(inv_matCh,obj.alpha,obj.beta);
            obj.compute_Chomog_Derivatives(physicalProblem.dim.nstre,physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus,x,interpolation,filter);
            gradient = obj.derivative_projection_Chomog(inv_matCh,obj.alpha,obj.beta,obj.Chomog_Derivatives,physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus,physicalProblem.dim.nstre);                                   
            if isempty(obj.h_C_0)
                obj.h_C_0 = costfunc;
            else
                costfunc = costfunc/abs(obj.h_C_0);
            end
            gradient = gradient/abs(obj.h_C_0);
            
            obj.value = costfunc;
            obj.gradient = gradient;
        end
    end
end