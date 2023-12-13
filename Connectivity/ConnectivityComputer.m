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
            d = obj.density.project('P0');
            
            %d.fValues = 
            
            d.fValues = 1-d.fValues; 
            d.fValues = round(d.fValues);              
          
            s.density = obj.density;
            s.mesh    = obj.mesh;
            s = StiffnessEigenModesComputer(s);
            [eigNeuman,eigDirichlet]  = s.compute();
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


            s.type = 'compliance';
            sh = ShapeFunctional.create(s);

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
            obj.levelSet.updateFunction();

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
            s.filterType = 'Lump';
            s.mesh  = obj.mesh;
            s.trial = P1Function.create(obj.mesh,1);
            f = Filter.create(s);
            dens = f.compute(obj.levelSet.fun,'QUADRATIC');
           %w    = max(0,min(1,1-dens));
           % w = 1 - dens;
           % w(:) = 1;
%            s.fValues = w;%floor(2*(w-0.5))+1;
%            s.mesh    = obj.mesh;
%            obj.density = P0Function(s);
            obj.density = dens;
        end





    

   
        
    end
    
end