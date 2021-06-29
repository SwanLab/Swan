classdef ConstitutiveTensorFromVademecum < VariableFromVademecum 
               
    
    properties (Access = protected)
       fieldName = 'Ctensor'; 
    end
    
    methods (Access = public)
        
        function [C,dC] = compute(obj,x)
            obj.computeParamsInfo(x);    
            obj.setValuesToInterpolator(x);
            [C,dC] = obj.computeValues();
        end
        
    end
    
    methods (Access = private)
        
        function [C,dC] = computeValues(obj)
            nstre = size(obj.values,1);            
            C  = zeros(nstre,nstre,obj.nPoints);
            dC = zeros(nstre,nstre,obj.nPoints,obj.nParams);
            for i = 1:nstre
                for j = 1:nstre
                    cv = squeeze(obj.values(i,j,:,:));                    
                    [c,dc] = obj.interpolator.interpolate(cv);
                    C(i,j,:) = c;
                    dC(i,j,:,1) = dc(:,1);
                    dC(i,j,:,2) = dc(:,2);
                end
            end            
        end
        
    end
    
end
