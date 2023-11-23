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
       materialInterpolation
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.filterCharacteristicFunction();

             obj.levelSet.getUnfittedMesh().plot()
             obj.density.plot()
            shading flat
            colormap('gray');
            colormap(flipud(gray));  
            colorbar
            figure

             obj.computeEigenValue()
          
            obj.computeCompliance();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end

        function computeEigenValue(obj)
             obj.createMaterialInterpolator();            
             obj.createBoundaryConditions();               
             obj.createStiffnessMatrix();
             obj.computeMassMatrix();
             [eigNeuman,eigDirichlet] = obj.obtainLowestEigenValues();
            
        end
        
        function computeCompliance(obj)
            obj.materialInterpolation = obj.computeMaterialInterpolation();
            s.mesh = obj.mesh;
            s.type    = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.createBoundaryConditionsForElasticity();
            s.interpolationType = 'LINEAR';
            fem = FEM.create(s);
            fem.solve();

        end

        function mat = createMaterial(obj)
            dens = 1- obj.density.fValues;
            mat  = obj.materialInterpolation.computeMatProp(dens);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mat.kappa;
            s.mu    = mat.mu;
            mat = Material.create(s);
            mat.compute(s);
        end
    
        function matInt = computeMaterialInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'SIMPALL';
            c.nElem = obj.mesh.nelem;
            c.dim = '2D';
            c.constitutiveProperties.rho_plus = 1;
            c.constitutiveProperties.rho_minus = 0;
            c.constitutiveProperties.E_plus = 1;
            c.constitutiveProperties.E_minus = 1e-3;
            c.constitutiveProperties.nu_plus = 1/3;
            c.constitutiveProperties.nu_minus = 1/3;

            matInt = MaterialInterpolation.create(c);
        end

  

        function bc = createBoundaryConditionsForElasticity(obj)
            bM = obj.mesh.createBoundaryMesh();
            
            dirichletBc.boundaryId=1;
            dirichletBc.dof=[1,2];
            dirichletBc.value=[0,0];
            newmanBc.boundaryId=2;
            newmanBc.dof=[2];
            newmanBc.value=[-1];

            [dirichlet,pointload] = obj.createBc(bM,dirichletBc,newmanBc);
            bc.dirichlet=dirichlet;
            bc.pointload=pointload;   
        end        

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
            pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end        


        function createMesh(obj)
            x1 = linspace(0,2,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh(s);            
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.ndim       = 2;
            s.fracRadius = 0.5;
            s.coord      = obj.mesh.coord;
            sD.type = 'LevelSet';
            sD.mesh = obj.mesh;
            sD.creatorSettings = s;
            sD.initialCase = 'circleInclusion';
            obj.levelSet   = DesignVariable.create(sD);


%             s.ndim       = 2;
%             s.widthH = 1;
%             s.widthV = 0.5;
%             s.coord      = obj.mesh.coord;
%             sD.type = 'LevelSet';
%             sD.mesh = obj.mesh;
%             sD.creatorSettings = s;
%             sD.initialCase = 'rectangleInclusion';
%             obj.levelSet   = DesignVariable.create(sD);


        end
 
        function filterCharacteristicFunction(obj)
            s.mesh = obj.mesh;
            s.quadratureOrder = [];
            s.femSettings.scale = 'MACRO';
            s.designVariable = obj.levelSet;
            %s.domainType = obj.mesh.type;
            f = Filter_PDE_LevelSet(s);
            dens = f.getP0fromP1([]);
           w    = max(0,min(1,1-dens));
           % w = 1 - dens;
           % w(:) = 1;
            s.fValues = w;%floor(2*(w-0.5))+1;
            s.mesh    = obj.mesh;
            obj.density = P0Function(s);
        end

        function createMaterialInterpolator(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPThermal';
            s.alpha0         = 1e-5;
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
            s.function = obj.density;
            s.quadratureOrder = 'QUADRATIC';
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            obj.Mmatrix = lhs.compute();       
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
        
    end
    
end