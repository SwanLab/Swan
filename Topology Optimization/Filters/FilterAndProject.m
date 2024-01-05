classdef FilterAndProject < handle

    properties (Access = private)
        filter
        projector
    end

    properties (Access = private)
        mesh
        type
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams);
            obj.createProjector(cParams);
        end

        function xF = compute(obj,fun)
            switch obj.type
                case 'DesignVariable'
                    xF = obj.computeDesignVariable(fun);
                case 'Gradient'
                    xF = obj.computeGradient(fun);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.type = cParams.type;
        end

        function createFilter(obj,s)
            s.filterType = 'PDE';
            obj.filter   = Filter.create(s);
        end

        function createProjector(obj,cParams)
            s             = cParams.femSettings;
            obj.projector = HeavisideProjector(s);
        end

        function xReg = computeDesignVariable(obj,fun)
            xFiltered = obj.filter.compute(fun);
            obj.projector.updateFilteredField(xFiltered);
            obj.projector.project();
            xReg = obj.projector.projectedField;
        end

        function xReg = computeGradient(obj,fun)
            x = fun.fValues;
            obj.projector.derive();
            dRhoBar_dRhoTilde = obj.projector.derivatedProjectedField;
            dC_dRhoTilde      = x.*dRhoBar_dRhoTilde;
            s.fValues         = dC_dRhoTilde;
            s.mesh            = obj.mesh;
            regFun            = P1Function(s);
            xReg              = obj.filter.getP1fromP0(regFun);
        end
    end
end