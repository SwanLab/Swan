classdef DesignVariable < handle & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        type
        nVariables                
    end
    
    properties (Access = public)
        value   
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
    
    methods (Access = public, Abstract)
        update(obj,value)
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end        
        
    end
    
    methods (Access = public)
                
        function restart(obj)
            obj.value = obj.valueOld;
            obj.alpha = obj.alphaOld;
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
            obj.type = cParams.type;
            obj.mesh = cParams.mesh; 
            cParams.scalarProductSettings.nVariables = obj.nVariables;
            obj.createScalarProduct(cParams);
        end
        
    end
    
    methods (Access = private)
        
        function createScalarProduct(obj,cParams)
            s = cParams.scalarProductSettings;
            s.femSettings.mesh = obj.mesh;
            obj.scalarProduct = ScalarProduct(s);        
        end
        
    end
    
end

