classdef StiffnessEigenModesComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        conductivity 
        boundaryConditions
        Kmatrix
        Mmatrix        
    end
    
    properties (Access = private)
        mesh
        density
    end
    
    methods (Access = public)
        
        function obj = StiffnessEigenModesComputer(cParams)
            obj.init(cParams)            
        end

        function [eigNeuman,eigDirichlet]  = compute(obj)
            obj.createMaterialInterpolator();
            obj.createBoundaryConditions();
            obj.createStiffnessMatrix();
            obj.computeMassMatrix();
            [eigNeuman,eigDirichlet] = obj.obtainLowestEigenValues();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.density = cParams.density;
            obj.mesh    = cParams.mesh;
        end


        function createMaterialInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.alpha0         = 1e-5;
            s.alpha1         = 1;
            s.density        = obj.density;
            a = MaterialInterpolation.create(s);
            obj.conductivity = a;            
%             sP.mesh          = obj.mesh;
%             sP.projectorType = 'P1D';
%             proj = Projector.create(sP);
%             a1   = proj.project(a);
        end    

        function createBoundaryConditions(obj)
            bMesh  = obj.mesh.createBoundaryMesh();
            allNodes = [];
            for i = 1:4
                nodes = bMesh{i}.globalConnec;
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

        function computeMassMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.function = obj.density;
            s.quadratureOrder = 'QUADRATIC';
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            obj.Mmatrix = lhs.compute();       
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
        
        
        function [eigLHSNewman,eigLHSDirichlet] = obtainLowestEigenValues(obj)
            K = obj.Kmatrix;
            M = obj.Mmatrix;
            bc  = obj.boundaryConditions;
            Kr = bc.fullToReducedMatrix(K);
            Mr = bc.fullToReducedMatrix(M);
            [V,eigLHSNewman]  = eigs(K,M,10,'smallestabs');
            [Vr,eigLHSDirichlet] = eigs(Kr,Mr,10,'smallestabs');

        %    [V,eigLHSNewman]  = eigs(K,[],10,'smallestabs');
         %   [Vr,eigLHSDirichlet] = eigs(Kr,[],10,'smallestabs');


            fV = zeros(size(V(:,1)));
            fV(obj.boundaryConditions.dirichlet,1) = obj.boundaryConditions.dirichlet_values;
            fV(obj.boundaryConditions.free,1) = Vr(:,1);
            s.fValues = fV;
            s.mesh    = obj.mesh;
            vV = P1Function(s);
            vV.plot()


            fV = V(:,2);
            s.fValues = fV;
            s.mesh    = obj.mesh;
            vN = P1Function(s);
            vN.plot()

        end   


        
    end
    
end