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

        function xReg = getP0fromP1(obj,x)
            xFiltered = obj.filter.getP0fromP1(x);
            obj.projector.updateFilteredField(xFiltered);
            obj.projector.project();
            xReg = obj.projector.projectedField;
        end

        function xReg = getP1fromP0(obj,x)
            obj.projector.derive();
            dRhoBar_dRhoTilde = obj.projector.derivatedProjectedField;
            dC_dRhoTilde      = x.*dRhoBar_dRhoTilde;
            xReg              = obj.filter.getP1fromP0(dC_dRhoTilde);
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

    end
    
end