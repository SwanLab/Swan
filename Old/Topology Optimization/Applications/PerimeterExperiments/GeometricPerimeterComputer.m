classdef GeometricPerimeterComputer < handle
    
   properties (Access = public)      
       perimeter      
   end
   
   properties (Access = private)
      unfittedMesh 
      circunferenceMesh
      integrator
      integrand
   end
    
   methods (Access = public)
      
       function obj = GeometricPerimeterComputer(cParams)
           obj.init(cParams);
       end
       
       function compute(obj)
            per = obj.unfittedMesh.computePerimeter;
            obj.perimeter = per;
       end
       
   end
   
   methods (Access = private)
       
       function init(obj,cParams)           
           obj.unfittedMesh = cParams.unfittedMesh;
       end
       
   end   
    
end