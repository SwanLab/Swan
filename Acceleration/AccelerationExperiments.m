classdef AccelerationExperiments < handle

    properties (Access = public)
        experiment
    end

    properties (Access = private)
        settings
        %
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = AccelerationExperiments(filename)
            obj.init(filename);
            obj.prepareOptimization();
            obj.createCaseAndCompute();
            obj.postprocessCase();
        end

    end

    methods (Access = private)

        function init(obj,filename)
            run(filename);
            obj.settings = s;
        end

        function prepareOptimization(obj)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createVolume();
            obj.createCost();
        end

        function createCaseAndCompute(obj)
            s.settings       = obj.settings;
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            obj.experiment   = AccelerationProblem(s);
            obj.experiment.createProblemAndCompute();
        end

        function createMesh(obj)
            [x1,x2] = obj.returnMeshCoordinates();
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);            
        end

        function [x1,x2] = returnMeshCoordinates(obj)
            switch obj.settings.geometryCase
                case {'CANTILEVER','ARCH','BRIDGE'}
                    x1 = linspace(0,2,50);
                    x2 = linspace(0,1.5,36);
                otherwise
                    error('Case not implemented yet')
            end
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(squeezeParticular(x(1,:,:),1)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);            
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;                        
            s.type    = 'Density';
            dens      = DesignVariable.create(s);   
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;            
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f.project('P1');
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createInterpolatedMaterial(f);
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.materialInterpolator        = obj.materialInterpolator;
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolume(obj)
            s.mesh     = obj.mesh;
            obj.volume = VolumeFunctional(s);
        end

        function createCost(obj)
            s.shapeFunctions = {obj.compliance,obj.volume};
            s.weights        = [1,obj.settings.lambda];
            obj.cost         = Cost(s);
        end

        function mat = createInterpolatedMaterial(obj,dens)
            mI   = obj.materialInterpolator;
            mat  = mI.computeConsitutiveTensor(dens);
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));

            sDir.direction = [1,2];
            sDir.value     = 0;
            sPL.direction  = 2;
            sPL.value      = -1;
           
            switch obj.settings.geometryCase
                case 'CANTILEVER'
                    isDir = @(coor) abs(coor(:,1)) == 0;
                    isPL  = @(coor) (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);
                case 'ARCH'
                    isDir = @(coor) abs(coor(:,2)) == 0 & (abs(coor(:,1))<=0.2*xMax | abs(coor(:,1))>=0.8*xMax);
                    isPL  = @(coor) abs(coor(:,2)) == 0 & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;
                case 'BRIDGE'
                    isDir = @(coor) abs(coor(:,2)) == 0;
                    isPL  = @(coor) abs(coor(:,2)) == yMax;
            end
            sDir.domain    = isDir;
            sPL.domain     = isPL;
            dirFun         = DirichletCondition(obj.mesh, sDir);
            plFun          = PointLoad(obj.mesh, sPL);
            s.dirichletFun = dirFun;
            s.pointloadFun = plFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function postprocessCase(obj)
            switch obj.settings.problemType
                case 'GENERAL'
                    obj.computeGeneralPlots();
                case 'TAU_BETA'
                otherwise 
            end
        end

        function computeGeneralPlots(obj)
            % legendNames = {'PG','PG-Polyak $\beta_{adapt.}$','PG-Polyak $\beta = 1$',...
            %     'PG-Nesterov $\beta_{adapt.}$  ','PG-Nesterov $\beta = 1$',...
            %     'MMA','MMA-Polyak $\beta_{adapt.}$','MMA-Polyak $\beta = 1$','Interpreter','Latex'};
            legendNames = {'PG','PG-Polyak $\beta_{adapt.}$','PG-Polyak $\beta = 1$',...
                'PG-Nesterov $\beta_{adapt.}$  ','PG-Nesterov $\beta = 1$',...
                'MMA','Interpreter','Latex'};
            minJ = inf;
            prob = obj.experiment.problem;
            for i = 1:numel(prob)
                minJ = min(minJ,min(prob{i}.J));
            end

            figure()
            hold on
            for i = 1:numel(prob)
                Jplot = prob{i}.J - minJ;
                hold on
                plot(Jplot,'-+','LineWidth',1.5)
            end
            xlabel('Iteration','interpreter','latex')
            ylabel('$\vert J - J^*\vert$','interpreter','latex')
            set(gca, 'YScale', 'log')
            box on
            grid minor
            legend(legendNames{:});
            set(gca,'FontSize',14,'TickLabelInterpreter','latex')
            hold off

            figure()
            hold on
            for i = 1:numel(prob)
                plot(prob{i}.incJ,'-+','LineWidth',1.5)
            end
            xlabel('Iteration','interpreter','latex')
            ylabel('$\Delta J$','interpreter','latex')
            set(gca, 'YScale', 'log')
            box on
            grid minor
            legend(legendNames{:});
            set(gca,'FontSize',14,'TickLabelInterpreter','latex')
            hold off

            if obj.settings.result
                m = obj.mesh;
                for i = 1:numel(prob)
                    TO = triangulation(m.connec,m.coord(:,1),m.coord(:,2),prob{i}.xFinal);
                    figure()
                    hold on
                    trimesh(TO,'facecolor','interp','edgecolor','none');
                    view(0,90)
                    colorbar
                    colormap(flipud(gray))
                    grid off
                    title(legendNames{i},'Interpreter','Latex')
                    box on
                    hold off
                end
            end
        end

    end

end