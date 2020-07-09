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
            obj.npnod = obj.mesh.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.Anodal2Gauss = obj.computeA();            
            cParams = SettingsMeshUnfitted(obj.domainType,obj.mesh);
            cParams.isInBoundary = false;
            obj.unfittedMesh = UnfittedMesh(cParams);
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
        
        function fInt = computeRHS(obj,ls,fNodes)
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                obj.unfittedMesh.compute(ls); 
                s.mesh = obj.unfittedMesh;
                s.type = 'Unfitted';
                int = Integrator.create(s);            
                fInt = int.integrateInDomain(fNodes);                    
            end
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