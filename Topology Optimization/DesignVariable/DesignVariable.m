classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        fun
        type        
    end

    properties (Access = protected)
        mesh
        nVariables
        isFixed
    end
    
    properties (Access = private)
        funOld
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end

    end

    methods (Access = public)
        
        function updateOld(obj)
            obj.funOld = obj.fun.copy();
        end

        function dX = computeIncrement(obj)
            dX = obj.fun-obj.funOld;
        end
        
        function norm = computeL2normIncrement(obj)
           incFun = obj.fun-obj.funOld;
           nIncX  = Norm.computeL2(obj.mesh,incFun);
           nX0    = Norm.computeL2(obj.mesh,obj.funOld);
           norm   = nIncX/nX0;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type = cParams.type;
            obj.mesh = cParams.mesh;
            obj.fun  = cParams.fun;
            if isfield(cParams,'isFixed')
              obj.isFixed = cParams.isFixed;
            end
        end
        
    end

end