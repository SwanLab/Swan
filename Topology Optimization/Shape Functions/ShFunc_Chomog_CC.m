classdef ShFunc_Chomog_CC < ShFunc_Chomog %%NOT WORKING%%
    properties (Access = private)
        Ch_star
        selectiveC_Cstar
    end
    methods
        function obj=ShFunc_Chomog_CC(settings)
            obj@ShFunc_Chomog(settings);
            obj.Ch_star=obj.compute_Ch_star(settings.TOL);
            obj.selectiveC_Cstar=settings.selectiveC_Cstar;
        end
        function computef(obj,x)
            obj.computePhysicalData(x);
            
            %Cost
            costfunc = obj.Chomog - obj.Ch_star;
            costfunc = obj.selectiveC_Cstar.*costfunc;
            
            %Gradient
            obj.compute_Chomog_Derivatives(x);
            DtC1 = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
            gradient = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
            for igaus=1:obj.physicalProblem.geometry.ngaus
                for a=1:obj.physicalProblem.dim.nstre
                    for b=1:obj.physicalProblem.dim.nstre
                        DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(a,b,igaus,:));
                        gradient(igaus,:) = gradient(igaus,:) + 2*costfunc(a,b)*DtC1(igaus,:);
                    end
                end
            end
            
            %Cost
            costfunc = sum(bsxfun(@times,costfunc(:),costfunc(:)));
            
            mass=obj.Msmooth;
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