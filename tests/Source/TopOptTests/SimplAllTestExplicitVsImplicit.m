classdef SimplAllTestExplicitVsImplicit < testShowingError
    
     properties (Access = protected)
       tol = 1e-12;        
     end    
     
     properties (Access = protected)
        rho
        matInterpExplicit
        matInterpImplicit
        materialPropertiesSettings        
        nElem
        nGaus
        dim
     end     
    
    methods (Access = protected)
        
        function computeError(obj)  
            mE = obj.matInterpExplicit;
            mI = obj.matInterpImplicit;      
            err(1) = norm(mE.mu(:) - mI.mu(:));
            err(2) = norm(mE.kappa(:) - mI.kappa(:));            
            err(3) = norm(mE.dmu(:) - mI.dmu(:));                        
            err(4) = norm(mE.dkappa(:) - mI.dkappa(:));                        
            obj.error = max(err);
        end
        
        function init(obj)
            obj.nElem = 800;
            obj.nGaus = 4;
            obj.rho = rand(obj.nElem,obj.nGaus);   
        end    
        
        function createProperties(obj)
            sC = SettingsConstitutiveProperties();            
            sC.E_plus   = rand(1);
            sC.E_minus  = rand(1);
            sC.nu_minus = rand(1);
            sC.nu_plus  = rand(1);                 
            obj.materialPropertiesSettings = sC;
        end
        
        function computeExplicitMaterialInterpolation(obj)            
            type = 'EXPLICIT';
            mI = obj.computeMaterialInterpolation(type);
            obj.matInterpExplicit = mI;
        end 
        
        function computeImplicitMaterialInterpolation(obj)            
            type = 'IMPLICIT';
            mI = obj.computeMaterialInterpolation(type);
            obj.matInterpImplicit = mI;
        end                  
        
        function p = computeMaterialInterpolation(obj,type)            
            sC = obj.materialPropertiesSettings;
            s = SettingsInterpolation();              
            s.constitutiveProperties = sC;            
            s.nElem = obj.nElem;   
            s.simpAllType = type;
            s.dim = obj.dim;
            mI = MaterialInterpolation.create(s);
            p = mI.computeMatProp(obj.rho);  
            obj.matInterpExplicit = p;
        end
        
    end   
    
end