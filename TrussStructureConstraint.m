classdef TrussStructureConstraint < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunction
    end
    
    methods (Access = public)
        
        function obj = TrussStructureConstraint(cParams)
            nCr = cParams.nConstraints;
            obj.value    = zeros(nCr,1);
            obj.gradient = zeros(length(cParams.designVariable),nCr);
            for iSF = nCr
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
            obj.value(iSF,1) = obj.shapeFunction{iSF}.value;
        end

        function computeGradient(obj,iSF)
            obj.shapeFunction{iSF}.computeGradient();
            obj.gradient(:,iSF) = obj.shapeFunction{iSF}.gradient;
        end
        
    end
    
    methods (Access = private)
        
        
    end

end
    