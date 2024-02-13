classdef HomogenizedMicrostructrueInterpolator < handle
    
    properties (Access = private)
        fileName
    end
    
    methods (Access = public)
        
        function obj = HomogenizedMicrostructrueInterpolator(cParams)
            obj.init(cParams)
        end

        function C = computeConsitutiveTensor(obj,x)
            C = obj.createMaterial(x);
        end

        function computeConstitutiveTensorDerivative(obj,x)
            
            for iVar = 1:length(x)
               dCi = dC(:,:,:,iVar);
               dCm{iVar} = obj.createMaterial(dCi);
            end            
        end
        
    end
 
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName = cParams.fileName;            
        end
      
        function m = createMaterial(obj,x)
            s.type        = 'ANISOTROPIC';
            s.fileName    = obj.fileName;            
            s.microParams = x;
            m = Material.create(s);   
        end
        
    end
    
    methods (Access = private, Static)
        
        function xR = reshapeDesignVariable(x)
            xV = x.value;
            nV = x.nVariables;
            nx = length(xV)/nV;
            xR = cell(nV,1);
            for ivar = 1:nV
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xR{ivar} = xV(i0:iF);
            end
        end
        
    end    
    
end