classdef FilterAndProject < Filter

    properties (Access = private)
        filter
        projector
    end
    
    methods (Access = public)

        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams);
            obj.createProjector(cParams);
        end

        function x_reg = getP0fromP1(obj,x)
            x_filtered = obj.computeFilter(x);
            obj.projector.updateFilteredField(x_filtered);
            obj.projector.project();
            x_reg = obj.projector.projectedField;
        end

        function x_reg = getP1fromP0(obj,x)
            obj.projector.derive();
            dRhoBar_dRhoTilde = obj.projector.derivatedProjectedField;
            dC_dRhoTilde = x.*dRhoBar_dRhoTilde;
            x_reg = obj.filter.getP1fromP0(dC_dRhoTilde);
        end
        
    end
    methods (Access = private)

        function createFilter(obj,s)
            s.filterType = 'PDE';
            obj.filter   = Filter.create(s);
        end

        function createProjector(obj,cParams)
            s             = cParams.femSettings;
            obj.projector = HeavisideProjector(s);
        end

        function x_filtered = computeFilter(obj,density)
            x_filtered = obj.filter.getP0fromP1(density);
        end

    end
end