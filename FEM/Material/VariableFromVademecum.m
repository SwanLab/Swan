classdef VariableFromVademecum < handle
    
    
    properties (Access = protected, Abstract)
        fieldName
    end
    
    properties (Access = protected)
        vadVariables
        interpolator
        values
        nPoints
        nParams
    end
    
    methods (Access = protected)
        
        function obj = VariableFromVademecum(cParams)
            obj.init(cParams);
            obj.obtainValues();
        end
        
        function setValuesToInterpolator(obj,x)
            obj.interpolator.setValues(x(:,1),x(:,2));
        end
        
        function computeParamsInfo(obj,x)
            [m,n] = size(x);
            obj.nPoints = m;
            obj.nParams = n;
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
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.vadVariables = cParams.vadVariables;
            obj.interpolator = cParams.interpolator;
        end
        
    end
       
end