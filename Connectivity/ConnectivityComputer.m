classdef ConnectivityComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       mesh
       levelSet
       density
       conductivity
       boundaryConditions
       Kmatrix
       Mmatrix
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.filterCharacteristicFunction();
            obj.createMaterialInterpolator();
            obj.createBoundaryConditions();            
            obj.createStiffnessMatrix();
            obj.computeMassMatrix();
            [eigLHS,eigLHSr] = obj.obtainLowestEigenValues();
            obj.density.plot()
            shading flat
            colormap('gray');
            colormap(flipud(gray));  
            colorbar
            figure
            obj.levelSet.getUnfittedMesh().plot()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function createMesh(obj)
            x1 = linspace(-1,1,40);
            x2 = linspace(-1,1,40);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh(s);            
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.ndim       = 2;
            s.fracRadius = 0.45;
            s.coord      = obj.mesh.coord;
            sD.type = 'LevelSet';
            sD.mesh = obj.mesh;
            sD.creatorSettings = s;
            sD.initialCase = 'circleInclusion';
            obj.levelSet   = DesignVariable.create(sD);
        end
 
        function filterCharacteristicFunction(obj)
            s.mesh = obj.mesh;
            s.quadratureOrder = [];
            s.femSettings.scale = 'MACRO';
            s.designVariable = obj.levelSet;
            %s.domainType = obj.mesh.type;
            f = Filter_PDE_LevelSet(s);
            dens = f.getP0fromP1([]);
            s.fValues = 1-dens;
            s.mesh    = obj.mesh;
            obj.density = P0Function(s);
        end

        function createMaterialInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.alpha0         = 1e-3;
            s.alpha1         = 1;
            s.density        = obj.density;
            a = MaterialInterpolation.create(s);
            obj.conductivity = a;
            
            sP.mesh          = obj.mesh;
            sP.projectorType = 'P1D';
            proj = Projector.create(sP);
            a1   = proj.project(a);
        end

        function createStiffnessMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';
            s.function        = obj.conductivity;
            s.type            = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            obj.Kmatrix = lhs.compute();
        end

        function computeMassMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';
            s.type            = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            obj.Mmatrix = lhs.compute();            
        end

        function [eigLHS,eigLHSr] = obtainLowestEigenValues(obj)
            K = obj.Kmatrix;
            M = obj.Mmatrix;
            bc  = obj.boundaryConditions;
            Kr = bc.fullToReducedMatrix(K);
            Mr = bc.fullToReducedMatrix(M);
            eigLHS  = eigs(K,M,10,'smallestabs');
            eigLHSr = eigs(Kr,Mr,10,'smallestabs');
        end

        function createBoundaryConditions(obj)
            bMesh  = obj.mesh.createBoundaryMesh();
            allNodes = [];
            for i = 1:4
                nodes = bMesh{1}.globalConnec;
                allNodes = [nodes;allNodes];
            end
            uNodes = unique(allNodes(:));

            bc{1}.ndimf = 1;
            bc{1}.ndofs = obj.mesh.nnodes;
            bc{1}.pointload = [];
            bc{1}.dirichlet(:,1) = uNodes;    
            bc{1}.dirichlet(:,2) = 1;    
            bc{1}.dirichlet(:,3) = 0;    

            s.scale   = 'MACRO';
            s.ndofs  = obj.mesh.nnodes;
            s.ndimf  = bc{1}.ndimf; 
            s.bc     = bc;
            bC = BoundaryConditions(s);
            obj.boundaryConditions = bC;
        end
        
    end
    
end