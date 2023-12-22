classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        type
        nVariables
        fun
    end
    
    properties (Access = public)
        alpha
        rho
    end
    
    properties (Access = private)
        valueOld
        alphaOld
    end
    
    properties (Access = protected)
       isFixed
    end
    
    methods (Access = public, Abstract)
        getVariablesToPlot(obj)
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function restart(obj)
            obj.update(obj.valueOld);
            obj.alpha = obj.alphaOld;
        end
        
        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            obj.fun   = P1Function(s);
        end
        
        function updateOld(obj)
            obj.valueOld = obj.fun.fValues;
            obj.alphaOld = obj.alpha;
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
        function norm = computeL2normIncrement(obj)
           l2Norm      = L2Norm(obj.mesh);
           x           = obj.fun.fValues;
           x0          = obj.valueOld;
           siF.fValues = x-x0;
           siF.mesh    = obj.mesh;
           incFun      = P1Function(siF);
           s0.fValues  = x0;
           s0.mesh     = obj.mesh;
           oldFun      = P1Function(s0);
           nIncX       = l2Norm.compute(incFun);
           nX0         = l2Norm.compute(oldFun);
           norm        = nIncX/nX0;
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