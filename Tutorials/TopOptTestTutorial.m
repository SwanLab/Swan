classdef TopOptTestTutorial < handle

    properties (Access = private)
        mesh
        filter
        designVariable
    end

    methods (Access = public)

        function obj = TopOptTestTutorial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();            
            obj.createFilter();
            obj.createElasticProblem();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);            
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);            
            s.fun     = aFun;
            s.mesh    = obj.mesh;                        
            s.type = 'Density';
            dens    = DesignVariable.create(s);   
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end        

        function createElasticProblem(obj)
            %Delete FEM
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            fem = ElasticProblem(s);
            fem.solve();
        end

       function mat = createMaterial(obj)
            matI = obj.computeMaterialInterpolation();           
            d = obj.designVariable.fun.project('P0');            
            dens    = d.fValues;
            mat     = matI.computeMatProp(dens);
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
            c.interpolation  = 'SIMPALL';
            c.nElem          = obj.mesh.nelem;
            c.dim            = '2D';
            
            cp.rho_plus = 1;
            cp.rho_minus = 0;
            cp.E_plus = 1;
            cp.E_minus = 1e-3;
            cp.nu_plus = 1/3;
            cp.nu_minus = 1/3;
            c.constitutiveProperties = cp;

            matInt = MaterialInterpolation.create(c);
        end       
        
        function bc = createBoundaryConditions(obj)
            bM = obj.mesh.createBoundaryMesh();

            dBC.boundaryId   = 1;
            dBC.dof          = [1,2];
            dBC.value        = [0,0];
            nBC.boundaryId   = 2;
            nBC.dof          = 2;
            nBC.value        = -1;

            [dirichlet,pointload] = obj.createBc(bM,dBC,nBC);
            bc.dirichlet=dirichlet;
            bc.pointload=pointload;
        end

       function [dirichlet,pointload] = createBc(obj,bMesh,dBC,nBC)
            dirichlet = obj.createBondaryCondition(bMesh,dBC);
            pointload = obj.createBondaryCondition(bMesh,nBC);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId = unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end        

    end

end