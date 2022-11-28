classdef TrussStructuresConstraint < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunction
    end
    
    methods (Access = public)
        
        function obj = TrussStructuresConstraint(cParams)
            for iSF = cParams.nConstraints
                cParams.type = cParams.constraintType{isF};
                obj.shapeFunction{iSF} = ShapeFunctional_Factory.create(cParams);
            end
        end

        function computeFunctionAndGradient(obj)
            for iSF = cParams.nConstraints
                obj.computeFunction(iSF);
                obj.computeGradient(iSF);
            end
        end

        function computeFunction(obj,iSF)
            obj.shapeFunction{iSF}.computeFunction();
            obj.value = obj.shapeFunction{iSF}.value;
        end

        function computeGradient(obj,iSF)
            obj.shapeFunction{iSF}.computeGradient();
            obj.gradient = obj.shapeFunction{iSF}.gradient;
        end
        
    end
    
    methods (Access = private)
        
        
    end

end
    