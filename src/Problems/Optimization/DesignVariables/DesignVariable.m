classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        fun
        type        
    end

    properties (Access = protected)
        nVariables
        isFixed
    end
    
    properties (Access = private)
        fValuesOld
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end

    end

    methods (Access = public)
        
        function recoverOld(obj)
            obj.fun.setFValues(obj.fValuesOld);
        end

        function updateOld(obj)
            obj.fValuesOld = obj.fun.fValues;
        end

        function norm = computeL2normIncrement(obj)
            funOld = obj.fun.copy();
            funOld.setFValues(obj.fValuesOld);
            incFun = obj.fun-funOld;
            nIncX  = Norm(incFun,'L2');
            nX0    = Norm(funOld,'L2');
            norm   = nIncX/nX0;
        end

    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type = cParams.type;
            obj.fun  = cParams.fun;
            if isfield(cParams,'isFixed')
              obj.isFixed = cParams.isFixed;
            end
        end
        
    end

end