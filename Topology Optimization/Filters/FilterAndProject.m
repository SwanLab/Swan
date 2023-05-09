classdef FilterAndProject < Filter
    properties (Access = public)
        filteredField
        projectedField
    end
    properties (Access = private)
        beta
        eta
        filter
        projector
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams)
            obj.defineProjectorSettings(cParams)
        end

        function compute(obj,density)
            obj.computeFilter(density);
            obj.createProjector()  
            obj.computeProjector()
        end
    end
    methods (Access = private)
        function createFilter(obj,cParams)
            s = cParams.filterParams;
            s.femSettings.mesh = s.mesh;
            s.designVariable = cParams.designVariable;
            obj.filter = Filter.create(s);
        end
        function defineProjectorSettings(obj,cParams)
            obj.eta = cParams.eta;
            obj.beta = cParams.beta;
        end 
        function computeFilter(obj,density)
            obj.filteredField = obj.filter.getP0fromP1(density);
        end
        function createProjector(obj)
            s.beta = obj.beta;
            s.eta = obj.eta;
            s.filteredField = obj.filteredField ;
            obj.projector = FieldProjector(s);
        end
        function computeProjector(obj)
            obj.projector.compute();
            obj.projectedField = obj.projector.projectedField;
        end 
    end
end