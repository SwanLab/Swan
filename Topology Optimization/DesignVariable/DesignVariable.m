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

        function x = obtainDomainFunction(designVariable)
            switch designVariable.type
                case 'Density'
                    x = designVariable.fun;
                case 'LevelSet'
                    x = designVariable.getCharacteristicFunction();
            end
        end

    end

    methods (Access = public)
        
        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            obj.fun   = P1Function(s);
        end
        
        function updateOld(obj)
            obj.funOld = obj.fun.copy();
        end
        
        function norm = computeL2normIncrement(obj)
           incFun = P1Function.substract(obj.fun,obj.funOld);
           nIncX  = Norm.computeL2(m,incFun);
           nX0    = Norm.computeL2(m,obj.funOld);
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