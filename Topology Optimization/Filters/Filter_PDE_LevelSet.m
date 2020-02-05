classdef Filter_PDE_LevelSet < Filter_PDE 
    
    properties(Access = private)
        integrator
        unfittedMesh
        domainType  
        interp
    end
    
    methods (Access = public)
        
        function obj = Filter_PDE_LevelSet(cParams)
            obj.init(cParams)
            obj.domainType = cParams.domainType;
            
            obj.diffReacProb.preProcess();
            obj.createQuadrature();
            obj.createInterpolation();            
            obj.computeGeometry();
            obj.nelem = obj.mesh.nelem;
            obj.npnod = obj.interp.npnod;           
            obj.ngaus = obj.quadrature.ngaus;
            obj.Anodal2Gauss = obj.computeA();
            
%             cParams = SettingsMeshUnfitted(obj.domainType,obj.mesh,obj.interpolation);
            cParams = SettingsMeshUnfitted(obj.domainType,obj.mesh);
            obj.unfittedMesh = UnfittedMesh(cParams);
            s.mesh = obj.unfittedMesh;
            s.type = obj.unfittedMesh.unfittedType;
            obj.integrator = Integrator.create(s);            
            obj.disableDelaunayWarning();                 
        end
        
        function preProcess(obj)
       
        end
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            F = ones(size(x));
            RHS = obj.computeRHS(x,F);
        end
        
        function RHS = integrate_function_along_facets(obj,x,F)
            RHS = obj.computeRHS(x,F);
        end
        
        function fInt = computeRHS(obj,levelSet,fNodes)
            obj.unfittedMesh.compute(levelSet);
            fInt = obj.integrator.integrate(fNodes);
        end        
        
    end
    
    methods (Access = private)
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.geometryType);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end        
          
        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');    
        end        
        
        function computeGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);        
        end
                
        
      function disableDelaunayWarning(obj)
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end        
          
    end
    
end