classdef FilterAndProject < Filter
    properties (Access = public)
        filteredField
        projectedField
    end
    properties (Access = private)
        eta
        beta
        filterIntermediate
        projector
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams)
            obj.defineProjectorSettings(cParams)
        end

        function x_reg = getP0fromP1(obj,x)
            obj.computeFilter(x);
            obj.createProjector(); 
            obj.projector.project();
            x_reg = obj.projector.projectedField;
        end

        function x_reg = getP1fromP0(obj,x)
            obj.projector.derive();
            dRhoBar_dRhoTilde = obj.projector.derivatedProjectedField;
            dC_dRhoTilde = x.*dRhoBar_dRhoTilde;
            x_reg = obj.filterIntermediate.getP1fromP0(dC_dRhoTilde);
        end

    end
    methods (Access = private)
        function createFilter(obj,s)
            s.filterType = 'P1';
            obj.filterIntermediate = Filter.create(s);
        end
        function defineProjectorSettings(obj,cParams)
            obj.eta  = cParams.femSettings.eta;
            obj.beta = cParams.femSettings.beta;
        end 
        function computeFilter(obj,density)
            obj.filteredField = obj.filterIntermediate.getP0fromP1(density);
        end
        function createProjector(obj)
            s.beta = obj.beta;
            s.eta = obj.eta;
            s.filteredField = obj.filteredField ;
            obj.projector = HeavisideProjector(s);
        end
    end
end