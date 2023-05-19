classdef FilterAndProject < Filter
    properties (Access = public)
        filteredField
        projectedField
    end
    properties (Access = private)
        eta
        beta
        filter
        projector
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams)
            obj.defineProjectorSettings(cParams)
        end

        function projectedField = getP0fromP1(obj,field)
            obj.computeFilter(field);
            obj.createProjector(); 
            obj.projector.project();
            projectedField = obj.projector.projectedField;
        end

        function dC_dRho = getP1fromP0(obj,dC_dRhoBar)
            obj.projector.deriveProjectedFieldRespectFilteredField();
            dRhoBar_dRhoTilde = obj.projector.derivatedProjectedField;
            dC_dRhoTilde = dC_dRhoBar.*dRhoBar_dRhoTilde;
            dC_dRho = obj.filter.getP1fromP0(dC_dRhoTilde);
        end

    end
    methods (Access = private)
        function createFilter(obj,s)
            s.filterType = 'P1';
            obj.filter = Filter.create(s);
        end
        function defineProjectorSettings(obj,cParams)
            obj.eta  = cParams.femSettings.eta;
            obj.beta = cParams.femSettings.beta;
        end 
        function computeFilter(obj,density)
            obj.filteredField = obj.filter.getP0fromP1(density);
        end
        function createProjector(obj)
            s.beta = obj.beta;
            s.eta = obj.eta;
            s.filteredField = obj.filteredField ;
            obj.projector = HeavisideProjector(s);
        end
    end
end