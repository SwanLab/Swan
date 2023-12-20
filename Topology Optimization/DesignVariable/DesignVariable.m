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
    
    properties (GetAccess = public, SetAccess = private)
        scalarProduct
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
           x  = obj.fun.fValues;
           x0 = obj.valueOld;
           incX  = x - x0;
           nIncX = obj.scalarProduct.computeSP_M(incX,incX);
           nX0   = obj.scalarProduct.computeSP_M(x0,x0);
           norm  = nIncX/nX0;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type    = cParams.type;
            obj.mesh    = cParams.mesh;
            if isfield(cParams,'isFixed')            
              obj.isFixed = cParams.isFixed;
            end
            obj.initValue(cParams);
            if isprop(cParams,'scalarProductSettings')  
                obj.createScalarProduct(cParams);
            end
        end
        
    end
    
    methods (Access = private)

        function initValue(obj,cParams)
            phiFun      = cParams.levelSetFunction;
            s.feFunType = class(phiFun);
            s.mesh      = obj.mesh;
            s.ndimf     = 1;
            obj.fun     = FeFunction.createEmpty(s);
            phi         = phiFun.fValues;
            switch obj.type
                case 'Density'
                    value = 1 - heaviside(phi);
                case 'LevelSet'
                    value = phi;
            end
            obj.fun.fValues = value;
        end
        
        function createScalarProduct(obj,cParams)
            s = cParams.scalarProductSettings;
            s.nVariables = obj.nVariables;
            s.femSettings.mesh = obj.mesh;
            obj.scalarProduct = ScalarProduct(s);
        end
        
    end 
end