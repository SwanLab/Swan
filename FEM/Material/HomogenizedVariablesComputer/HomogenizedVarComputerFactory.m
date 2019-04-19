classdef HomogenizedVarComputerFactory < handle
    
   methods (Access = public, Static)
       
       function h = create(cParams)
           switch cParams.type
               case 'ByVademecum'
                   s = SettingsHomogenizedVarComputerFromVademecum();
                   s.fileName = cParams.vademecumFileName;
                   h = HomogenizedVarComputerFromVademecum(s);                   
               case 'ByInterpolation'
                   s = SettingsHomogenizedVarComputerFromInterpolation();
                   s.interpolation          = cParams.interpolation;
                   s.dim                    = cParams.dim;
                   s.typeOfMaterial         = cParams.typeOfMaterial;
                   s.constitutiveProperties = cParams.constitutiveProperties;  
                   s.nelem                  = cParams.nelem;
                   h = HomogenizedVarComputerFromInterpolation(s);
           end
           
       end       
       
   end
    
end