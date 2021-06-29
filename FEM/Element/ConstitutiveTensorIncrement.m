classdef ConstitutiveTensorIncrement < handle
    
    properties
        Cmat
        changedElements
        TOL = 1e-15
        dim
    end
    
    
    
    methods
        function obj = ConstitutiveTensorIncrement(dim)
            obj.dim = dim;
        end
        
        function  obtainChangedElements(obj,CmatNew)
            if ~isempty(obj.Cmat)
                
                dC = obj.Cmat - CmatNew;
                
                dCnorm = squeeze(max(max(abs(dC),[],1),[],2));
                CmatNewNorm = max(abs(CmatNew(:)));
                
                obj.changedElements.values  = (dCnorm)/CmatNewNorm > obj.TOL;
                obj.changedElements.number  = sum(obj.changedElements.values);
                obj.changedElements.Ratio   = obj.changedElements.number/size(obj.changedElements.values,1);
                obj.changedElements.Entries = repmat(obj.changedElements.values,obj.dim.ndofPerElement*obj.dim.ndofPerElement,1);
                
            else
                obj.changedElements.values  = true(obj.dim.nelem,1);
                obj.changedElements.number  = sum(obj.changedElements.values);
                obj.changedElements.Ratio   = obj.changedElements.number/size(obj.changedElements.values,1);
                
                %xobj.changedElements.Ratio   = 1;
                obj.changedElements.Entries = repmat(obj.changedElements.values,obj.dim.ndofPerElement*obj.dim.ndofPerElement,1);
            end
            
            obj.Cmat = CmatNew;
            
        end
        
        
    end
    
    
    
    
    
end

