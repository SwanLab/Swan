classdef DesignVariable < handle
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        type
        nVariables        
        value           
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
            obj.value = value;
            if ~isempty(obj.isFixed)
               obj.value(obj.isFixed.nodes) = obj.isFixed.values;
            end
        end        
        
        function updateOld(obj)
            obj.valueOld = obj.value;
            obj.alphaOld = obj.alpha;
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
        function norm = computeL2normIncrement(obj)
           x = obj.value;
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
            obj.isFixed = cParams.isFixed;
            obj.initValue(cParams);
            obj.createScalarProduct(cParams);
        end
        
    end
    
    methods (Access = private)
        
        function initValue(obj,cParams)
            if isempty(cParams.value)
                obj.value = ones(size(obj.mesh.coord,1),1);
            else
                obj.value = cParams.value;
            end
        end
        
        function createScalarProduct(obj,cParams)
            s = cParams.scalarProductSettings;
            s.nVariables = obj.nVariables;
            s.femSettings.mesh = obj.mesh;
            obj.scalarProduct = ScalarProduct(s);        
        end
        
    end

    
end

