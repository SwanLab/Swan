classdef GeometricPerimeterComputer < handle
    
   properties (Access = public)      
       perimeter      
   end
   
   properties (Access = private)
      designVariable 
      circunferenceMesh
      integrator
      integrand
   end
    
   methods (Access = public)
      
       function obj = GeometricPerimeterComputer(cParams)
           obj.init(cParams);
       end
       
       function compute(obj)
            obj.createCircunferenceMesh();           
            obj.createIntegrator();
            obj.integrateOneFunctionInCircunferenceMesh();          
            obj.perimeter = sum(obj.integrand);           
       end
       
   end
   
   methods (Access = private)
       
       function init(obj,cParams)           
           obj.designVariable = cParams.designVariable;
       end
       
       function createCircunferenceMesh(obj)
            s = obj.createMeshUnfittedSettings();
            cMesh = UnfittedMesh(s);
            cMesh.compute(obj.designVariable.value);           
            obj.circunferenceMesh = cMesh;
       end
       
       function s = createMeshUnfittedSettings(obj)
            mBackground   = obj.designVariable.mesh;
            sM.unfittedType            = 'BOUNDARY';
            sM.backgroundMesh = mBackground.innerMeshOLD;
            sM.boundaryMesh   = mBackground.boxFaceMeshes;
            s = SettingsMeshUnfitted(sM);           
       end
       
       function createIntegrator(obj)
            s.mesh  = obj.circunferenceMesh;
            s.type  = 'COMPOSITE';
            sc.mesh = obj.circunferenceMesh.innerCutMesh;
            sc.type = 'CutMesh';
            s.compositeParams{1} = sc;
            obj.integrator = Integrator.create(s);
       end
       
       function integrateOneFunctionInCircunferenceMesh(obj)
            npnod = obj.circunferenceMesh.meshBackground.npnod;
            f = ones(npnod,1);
            obj.integrand = obj.integrator.integrateAndSum(f);           
       end
       
   end   
    
end