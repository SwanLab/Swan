classdef ElasticProblemDavid < handle
    
    properties (Access = public)
        uFun
        sigmaFun
        strainFun
        mesh
        young 
        poisson
        material
        bc

        stiffness
        LHS
        RHS
    end
    
    properties (Access = private)
    end

    methods (Access = public)
        
        function obj = ElasticProblemDavid()
            obj.init()
            obj.createMesh()
            obj.computeElasticProperties()
            obj.createMaterial()
            obj.createBoundaryConditions()
            obj.displacementFunction()
            obj.createStiffness()
            obj.createForce()
            obj.computeNewDisplacement()
            obj.computeStressAndStrain()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end

        function createMesh(obj)
            %Creo malla 2x1x1 
            obj.mesh = HexaMesh(2,1,1,2,1,1);
        end

        function computeElasticProperties(obj)
            poisson = 0.33;
            young = 70;
            obj.poisson = ConstantFunction.create(poisson,obj.mesh);
            obj.young = ConstantFunction.create(young,obj.mesh);
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young = obj.young;
            s.poisson = obj.poisson;
            obj.material = Material.create(s);
        end

        function createBoundaryConditions(obj)
            
            %AIXO HO HE COPIAT PARCIALMENT
            right    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMax    = max(obj.mesh.coord(:,3));
             
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight = @(coor)  abs(coor(:,1))==right;


            sDir{1}.domain    = @(coor) isLeft(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isRight(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;
            
                    sPL{2}.domain    = @(coor) isRight(coor);
                    sPL{2}.direction = 3;
                    sPL{2}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;


            obj.bc = BoundaryConditions(s);



        end

        function displacementFunction(obj)
            nDimF = obj.mesh.ndim; 
            obj.uFun = LagrangianFunction.create(obj.mesh,nDimF,'P1'); 
        end

        function createStiffness(obj)

            C = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            obj.stiffness = IntegrateLHS (f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            
        end

        function createForce(obj)
            rhs = zeros(obj.uFun.nDofs,1);

            for i = 1:numel(obj.bc.tractionFun)
                rhs = rhs + obj.bc.tractionFun(i).computeRHS(obj.uFun); % això és f
            end

            % SolverType --> Reduced
                % rhs = f-Klr*ur
                R = -obj.stiffness(:,obj.bc.dirichlet_dofs)*obj.bc.dirichlet_vals; % això és Klr*ur
                obj.RHS = rhs + R;
        end

        function computeNewDisplacement(obj)
            dofs = 1:size(obj.stiffness,1);
            free_dofs = setdiff(dofs, obj.bc.dirichlet_dofs);

            obj.LHS = obj.stiffness(free_dofs, free_dofs);
            obj.RHS = obj.RHS(free_dofs); %obtingut a createForce

            u(free_dofs) = obj.LHS\obj.RHS; 
            u(obj.bc.dirichlet_dofs) = obj.bc.dirichlet_vals;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.uFun.setFValues(uSplit);
            
        end

        function obj = computeStressAndStrain(obj)
            obj.strainFun = SymGrad(obj.uFun);
            C = obj.material;
            obj.sigmaFun = DDP(C,obj.strainFun); 
            


        end
        
    end
    
end