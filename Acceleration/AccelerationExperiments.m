classdef AccelerationExperiments < handle

    properties (Access = public)
        experiment
    end

    properties (Access = public)
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
        solverTol
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
            obj.createSolverTolerance();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createVolume();
            obj.createCost();
        end

        function createCaseAndCompute(obj)
            s.settings       = obj.settings;
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.solverTol      = obj.solverTol;
            obj.experiment   = AccelerationProblem(s);
            tInit = tic;
            obj.experiment.createProblemAndCompute();
            toc(tInit)
        end

        function createMesh(obj)
            switch obj.settings.geometryType
                case '2D'
                    [x1,x2] = obj.returnMeshCoordinates();
                    [xv,yv] = meshgrid(x1,x2);
                    [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
                    s.coord  = V(:,1:2);
                    s.connec = F;
                    obj.mesh = Mesh.create(s);
                case '3D'
                    [x1,x2,x3] = obj.returnMesh3DCoordinates();
                    d1 = 60;%40;%
                    d2 = 30;%20;%
                    d3 = 30;%20;%
                    obj.mesh = TetraMesh(x1,x2,x3,d1,d2,d3);
            end
                        
        end

        function createSolverTolerance(obj)
            s.solver      = obj.settings.solverType;
            s.tolMax      = obj.settings.maxTol;
            s.tolMin      = obj.settings.minTol;
            obj.solverTol = ConjugateGradientToleranceCalculator(s);
        end

        function [x1,x2] = returnMeshCoordinates(obj)
            switch obj.settings.geometryCase
                case {'CANTILEVER','ARCH','BRIDGE'}
                    x1 = linspace(0,2,80);
                    x2 = linspace(0,1.5,60);
                    % x1 = linspace(0,2,40);
                    % x2 = linspace(0,1.5,20);
                otherwise
                    error('Case not implemented yet')
            end
        end

        function [x1,x2,x3] = returnMesh3DCoordinates(obj)
            switch obj.settings.geometryCase
                case {'CANTILEVER3D','ARCH3D','CANTILEVER_2'}
                    x1 = 2;
                    x2 = 0.75;
                    x3 = 1.5;
                case 'BRIDGE3D'
                    x1 = 2;
                    x2 = 1.2;
                    x3 = 1.5;
            end
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
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

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f.project('P1');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = obj.settings.geometryType;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            s.material = obj.createMaterial();
            s.dim      = obj.settings.geometryType;
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = obj.settings.solverType;
            s.matrixFree         = obj.settings.matrixFree;
            s.solverTol          = obj.solverTol;
            p.maxIters = 5e3;
            % if obj.settings.matrixFree p.maxIters = 5e3; 
            % else p.maxIters = 5e3; end       
            p.displayInfo        = true;
            s.solverParams       = p;
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstitutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolume(obj)
            s.mesh     = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.volume = VolumeFunctional(s);
        end

        function createCost(obj)
            s.shapeFunctions = {obj.compliance,obj.volume};
            s.weights        = [1,obj.settings.lambda];
            s.Msmooth        = obj.createMassMatrix();
            obj.cost         = Cost(s);
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;     
        end

        function mat = createInterpolatedMaterial(obj,dens)
            mI   = obj.materialInterpolator;
            mat  = mI.computeConsitutiveTensor(dens);
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
           
            switch obj.settings.geometryCase
                case 'CANTILEVER'
                    isDir = @(coor) abs(coor(:,1)) == 0;
                    isPL  = @(coor) (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    sPL.direction  = 2;
                    sPL.value      = -1;
                case 'ARCH'
                    isDir = @(coor) abs(coor(:,2)) == 0 & (abs(coor(:,1))<=0.2*xMax | abs(coor(:,1))>=0.8*xMax);
                    isPL  = @(coor) abs(coor(:,2)) == 0 & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    sPL.direction  = 2;
                    sPL.value      = -1;
                case 'BRIDGE'
                    isDir = @(coor) abs(coor(:,2)) == 0;
                    isPL  = @(coor) abs(coor(:,2)) == yMax;
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    sPL.direction  = 2;
                    sPL.value      = -1;
                case 'CANTILEVER3D'
                    zMax  = max(obj.mesh.coord(:,3));
                    isDir = @(coor)  abs(coor(:,1))==0;
                    isPL  = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.45*yMax & (abs(coor(:,2))<=0.55*yMax) ...
                        & abs(coor(:,3))>=0.45*zMax & abs(coor(:,3))<=0.55*zMax);
                    sDir.direction = [1,2,3];
                    sDir.value     = 0;
                    sPL.direction  = 3;
                    sPL.value      = -1;
                case 'CANTILEVER_2'
                    isDir = @(coor)  abs(coor(:,1))==0;
                    isPL  = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,3))==0);
                    sDir.direction = [1,2,3];
                    sDir.value     = 0;
                    sPL.direction  = 3;
                    sPL.value      = -1;
                case 'BRIDGE3D'
                    zMax  = max(obj.mesh.coord(:,3));
                    isDir = @(coor)  abs(coor(:,3))==0;
                    isPL  = @(coor)  abs(coor(:,3))==zMax & abs(coor(:,2))>=0.15*yMax & abs(coor(:,2))<=0.85*yMax;
                    sDir.direction = [1,2,3];
                    sDir.value     = 0;
                    sPL.direction  = 3;
                    sPL.value      = -1;
                case 'ARCH3D'
                    isDir = @(coor) abs(coor(:,3)) == 0 & (abs(coor(:,1))<=0.2*xMax | abs(coor(:,1))>=0.8*xMax);
                    isPL  = @(coor) abs(coor(:,3)) == 0 & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;
                    sDir.direction = [1,2,3];
                    sDir.value     = 0;
                    sPL.direction  = 3;
                    sPL.value      = -1;
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
                case 'SOLVER'
                    obj.computeSolverPlots();
                case 'GENERAL'
                    obj.computeGeneralPlots();
                case 'TAU_BETA'
                    obj.computeTauBetaPlots();
                otherwise 
            end
        end

        function computeSolverPlots(obj)
            p = obj.experiment.problem;
            if obj.settings.solverType == "CONJUGATE GRADIENT"
                nTols = obj.physicalProblem.getSolverTols();
                nIters = obj.physicalProblem.getSolverIters();

                t  = nTols(1:2:end);
                it = nIters(1:2:end);
                disp('Total CG iterations = ' + string(sum(it)));

                figure()
                plot(1:numel(t),t,'k')
                xlabel('TO iteration','interpreter','latex')
                ylabel('CG tolerance $\epsilon$','interpreter','latex')
                set(gca,'FontSize',14,'TickLabelInterpreter','latex')
                set(gca, 'YScale', 'log')
                grid minor
                xlim([1 numel(t)])
                ylim([min(t) max(t)*1.2])
                box on

                figure()
                bar(1:numel(it),it)
                xlabel('TO iteration','interpreter','latex')
                ylabel('CG iterations to converge','interpreter','latex')
                set(gca,'FontSize',14,'TickLabelInterpreter','latex')
                xlim([1 numel(it)])
                box on
            end

            figure()
            plot(p.costFields(1,:),'k')
            xlabel('TO iteration','interpreter','latex')
            ylabel('Compliance','interpreter','latex')
            set(gca,'FontSize',14,'TickLabelInterpreter','latex')
            grid minor
            xlim([1 numel(p.costFields(1,:))])
            box on

            figure()
            plot(p.costFields(2,:),'k')
            xlabel('TO iteration','interpreter','latex')
            ylabel('Volume','interpreter','latex')
            set(gca,'FontSize',14,'TickLabelInterpreter','latex')
            grid minor
            xlim([1 numel(p.costFields(1,:))])
            box on
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
                    axis equal
                    hold off
                end
            end

            if obj.settings.saveResults
                d = obj.designVariable;
                for i = 1:numel(prob)
                    d.update(prob{i}.xFinal);
                    d.fun.print('x')
                end
            end
        end

        function computeTauBetaPlots(obj)
            st = obj.settings;
            p = obj.experiment;
            [betaGrid,tauGrid] = meshgrid(st.beta,st.tau);
            bestOnes = zeros(size(p.problem,1),3);
            for i = 1:size(p.problem,1)
                [val,pos]     = min(p.problem(i,:));
                bestOnes(i,:) = [betaGrid(1,pos),tauGrid(i,1),val+1e3];
            end

            figure()
            hold on
            s = surf(betaGrid,tauGrid,p.problem);
            xlabel('Momentum term ($\beta$)','Interpreter','latex')
            ylabel('Line search ($\tau$)','Interpreter','latex')
            s.EdgeColor = "none";
            s.FaceColor = "interp";
            c = colorbar;
            colormap(flipud(jet));
            c.TickLabelInterpreter = 'latex';
            clim([min(min(p.problem)),st.maxIter]);
            view(0,90)
            h = plot3(bestOnes(:,1),bestOnes(:,2),bestOnes(:,3),'o','Color','black');
            set(h, 'MarkerFaceColor', 'k');
            set(gca,"TickLabelInterpreter",'latex','FontSize',14)
            xlim([st.beta(1) st.beta(end)])
            ylim([st.tau(1) st.tau(end)])
            box on
            hold off
            
        end

    end

end