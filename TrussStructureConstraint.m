classdef TrussStructureConstraint < handle
    
    properties (Access = public)
        value
        gradient
        nSF
    end

    properties (Access = private)
        shapeFunction
        nConstraints
    end
    
    methods (Access = public)
        
        function obj = TrussStructureConstraint(cParams)
            nCr = cParams.nConstraints;
            obj.nConstraints = nCr;
            obj.value    = zeros(nCr,1);
            obj.gradient = zeros(length(cParams.designVariable),nCr);
            for iSF = 1:nCr
                cParams.type = cParams.ctype{iSF};
                f = ShapeFunctional_Factory();
                obj.shapeFunction{iSF} = f.create(cParams);
            end
        end

        function computeFunctionAndGradient(obj)
            for iSF = obj.nConstraints
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
    