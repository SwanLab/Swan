classdef HomogenizedVarComputerFromInterpolation ...
         < HomogenizedVarComputer
    
    properties (Access = private)
        interpolation
        material
    end
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromInterpolation(cParams)
            obj.createMaterialInterpolation(cParams);
            s.nelem = cParams.nelem;
            s.ptype = 'ELASTIC';
            s.pdim  = cParams.dim;
            obj.material = Material.create(s);
        end
        
        function computeMatProp(obj,rho)
           mProps = obj.interpolation.computeMatProp(rho);
           obj.material.setProps(mProps);
           obj.C  = obj.material.C;
           obj.dC = mProps.dC;           
        end        
        
    end
    
    methods (Access = private)
        
        function createMaterialInterpolation(obj,cParams)
            s = SettingsInterpolation();
            s.interpolation          = cParams.interpolation;
            s.dim                    = cParams.dim;
            s.typeOfMaterial         = cParams.typeOfMaterial;
            s.constitutiveProperties = cParams.constitutiveProperties;            
            int = Material_Interpolation.create(s);        
            obj.interpolation = int;            
        end
        
    end
    
end