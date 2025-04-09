classdef StiffnessEigenModesComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        conductivity 
        massInterpolator
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
            obj.createBoundaryConditions();    
            obj.createConductivityInterpolator();   
            obj.createMassInterpolator();                        
        end

        function [lambda,dlambda]  = computeFunctionAndGradient(obj,dens)
            obj.density = dens;
            alpha  = obj.conductivity.fun;
            dalpha = obj.conductivity.dfun;
            m      = obj.massInterpolator.fun;
            dm     = obj.massInterpolator.dfun;            
            K  = obj.createStiffnessMatrixWithFunction(alpha);
            M  = obj.computeMassMatrixWithFunction(m);
            dK = obj.createStiffnessMatrixWithFunction(dalpha);
            dM = obj.computeMassMatrixWithFunction(dm);            
            [lambda,phi] = obj.obtainLowestEigenValuesAndFunction(K,M);
            dlambda = phi'*dK*phi - lambda*phi'*dM*phi;
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
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
        
        function createConductivityInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.f0   = 1e-5;
            s.f1   = 1;
            s.pExp = 8;
            a = MaterialInterpolation.create(s);
            obj.conductivity = a;            
        end            

        function createMassInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.f0   = 1e-5;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolation.create(s);
            obj.massInterpolator = a;            
        end            

        function K = createStiffnessMatrixWithFunction(obj,fun)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';
            s.function        = obj.createDomainFunction(fun);
            s.type            = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            K = obj.boundaryConditions.fullToReducedMatrix(K);
        end

        function f = createDomainFunction(obj,fun)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xV);
            s.mesh      = obj.mesh;
            f = DomainFunction(s);
        end

        function fV = createConductivityAsDomainFunction(obj,fun,xV)
            densV = obj.density.evaluate(xV);
            fV = fun(densV);
        end


        function M = computeMassMatrixWithFunction(obj,fun)
            s.test  = P1Function.create(obj.mesh,1); 
            s.trial = P1Function.create(obj.mesh,1); 
            s.mesh  = obj.mesh;
            s.function = obj.density;
            s.quadratureOrder = 'QUADRATIC';
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();   
            M = obj.boundaryConditions.fullToReducedMatrix(M);
        end       
                
        function [eigV1,eigF1] = obtainLowestEigenValuesAndFunction(obj,K,M)
            [eigF,eigV] = eigs(K,M,2,'smallestabs');
            eigV1 = eigV(1);
            eigF1 = eigF(1);
            

        %    [V,eigLHSNewman]  = eigs(K,[],10,'smallestabs');
         %   [Vr,eigLHSDirichlet] = eigs(Kr,[],10,'smallestabs');


            % fV = zeros(size(V(:,1)));
            % fV(obj.boundaryConditions.dirichlet,1) = obj.boundaryConditions.dirichlet_values;
            % fV(obj.boundaryConditions.free,1) = Vr(:,1);
            % s.fValues = fV;
            % s.mesh    = obj.mesh;
            % vV = P1Function(s);
            % vV.plot()
            % 
            % 
            % fV = V(:,2);
            % s.fValues = fV;
            % s.mesh    = obj.mesh;
            % vN = P1Function(s);
            % vN.plot()
    
        end   


        
    end
    
end