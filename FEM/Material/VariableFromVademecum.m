classdef VariableFromVademecum < handle
    
    
    properties (Access = protected, Abstract)
        fieldName
    end
    
    properties (Access = protected)
        vadVariables
        interpolator
        values
        np
        indexMx
        indexMy
    end    
    
    methods (Access = protected)
        
        function obj = VariableFromVademecum(cParams)
            obj.init(cParams);                      
            obj.obtainValues();
        end        
        
        function init(obj,cParams)
            obj.vadVariables = cParams.vadVariables;
            obj.interpolator = cParams.interpolator;
        end
        
        function obtainValues(obj)
            var = obj.vadVariables.variables;
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.mxV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    v(:,:,imx,imy) = var{imx,imy}.(obj.fieldName);
                end
            end
            obj.values = v;
        end
        
        function computeParamsInfo(obj,x)
            obj.np = length(x)/2;
            obj.indexMx = 1:obj.np;
            obj.indexMy = obj.np+1:2*obj.np;
        end
        
        function setValuesToInterpolator(obj,x)
            obj.interpolator.setValues(x(obj.indexMx,1),x(obj.indexMy,1));               
        end        
        
        
    end
    
end