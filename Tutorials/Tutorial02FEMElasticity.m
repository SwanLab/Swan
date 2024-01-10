classdef Tutorial02FEMElasticity < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
    end

    methods (Access = public)

        function obj = Tutorial02FEMElasticity()
            obj.init();
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
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

        function computeElasticProperties(obj)
            E0  = 1e-3;
            nu0 = 1/3;
            E   = AnalyticalFunction.create(@(x) E0*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu  = AnalyticalFunction.create(@(x)nu0*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            obj.young = E;
            obj.poisson = nu;
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.bc = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
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