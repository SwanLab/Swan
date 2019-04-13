classdef ConstitutiveTensorFromVademecum < VariableFromVademecum
    
    properties (Access = protected)
       fieldName = 'Ctensor'; 
    end
    
    methods (Access = public)
        
        function [C,dC] = computeCtensor(obj,x)
            obj.computeParamsInfo(x);    
            obj.setValuesToInterpolator(x);
            
            nstre = size(obj.values,1);            
            C  = zeros(nstre,nstre,obj.np);
            dC = zeros(nstre,nstre,2*obj.np);
            for i = 1:nstre
                for j = 1:nstre
                    cv = squeeze(obj.values(i,j,:,:));                    
                    [c,dc] = obj.interpolator.interpolate(cv);
                    C(i,j,:) = c;
                    dC(i,j,obj.indexMx) = dc(:,1);
                    dC(i,j,obj.indexMy) = dc(:,2);
                end
            end
        end
        
    end
    
end
