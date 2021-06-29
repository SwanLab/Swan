classdef Filter_PDE_LevelSet < Filter_PDE
    
    properties (Access = public)
        unfittedMesh        
    end
    
    properties(Access = private)
        integrator
        domainType
        interp
        levelSet
    end
    
    methods (Access = public)
        
        function obj = Filter_PDE_LevelSet(cParams)
            obj.init(cParams)
            obj.levelSet = cParams.designVariable;
            obj.epsilon = cParams.mesh.computeMeanCellSize();            
            obj.domainType = cParams.domainType;            
            obj.diffReacProb.preProcess();
            obj.createQuadrature();
            obj.createInterpolation();
            obj.computeGeometry();
            obj.nelem = obj.mesh.nelem;
            obj.npnod = obj.mesh.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.Anodal2Gauss = obj.computeA(); 
        end
        
        function preProcess(obj)
            preProcess@Filter(obj)                                    
            obj.Anodal2Gauss = obj.computeA();  
            obj.diffReacProb.setEpsilon(obj.epsilon);
            obj.computeLHS();
        end
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            ls = obj.levelSet.value;
            F = ones(size(ls));
            RHS = obj.computeRHS(F);
        end
        
        function RHS = integrate_function_along_facets(obj,F)
            RHS = obj.computeRHSinBoundary(F);
        end
        
    end
    
    methods (Access = private)
        
        function fInt = computeRHS(obj,fNodes)
            ls = obj.levelSet.value;
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                uMesh = obj.levelSet.getUnfittedMesh();
                s.mesh = uMesh;
                s.type = 'Unfitted';
                int = Integrator.create(s);            
                fInt = int.integrateInDomain(fNodes);                    
            end
        end    
        
        function fInt = computeRHSinBoundary(obj,fNodes)
            ls = obj.levelSet.value;
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                uMesh = obj.levelSet.getUnfittedMesh();
                s.mesh = uMesh;
                s.type = 'Unfitted';
                int = Integrator.create(s);            
                fInt = int.integrateInBoundary(fNodes);                    
            end
        end             
        
        
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
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
        
    end
    
end