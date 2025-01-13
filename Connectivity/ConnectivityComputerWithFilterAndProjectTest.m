classdef ConnectivityComputerProjectionTest < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        levelSet
        characteristicFunction
        designVariable
        materialInterpolation
        filter
        filterConnect
        filterAdjointConnect
        monitoringEigenModes
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputerProjectionTest()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable(); 
%             obj.createAnalyticalEigenModes(0,0);
%             obj.createAnalyticalEigenModes(1,0);
%             obj.createAnalyticalEigenModes(0,1);
%             obj.createAnalyticalEigenModes(1,1);
            figure;
%             analytical = [0.0, pi^2*(1/0.5^2), pi^2*(1/0.5^2), pi^2*(2/0.5^2), pi^2*(4/0.5^2), pi^2*(4/0.5^2), pi^2*(5/0.5^2), pi^2*(5/0.5^2), pi^2*(8/0.5^2), pi^2*(9/0.5^2), pi^2*(9/0.5^2), pi^2*(10/0.5^2), pi^2*(10/0.5^2), pi^2*(13/0.5^2), pi^2*(13/0.5^2), pi^2*(18/0.5^2)];
%             analytical = [0.0, 0.0, pi^2, 4*pi^2, 9*pi^2, 100*pi^2/9];
            analytical = [0.0, 4*pi^2,  4*pi^2,  8*pi^2,  16*pi^2,  16*pi^2 ];
            n = 6;
            betas = 1.0:5:300;
            j = 1;
            eigF1 = []; 
            eigF2 = [];
            eigF3 = [];
            eigF4 = [];
            for type = ["FP", "P"]
                for eta = [0.0, 0.5, 1.0]
                    eigVs= [];
                    eigFs= [];   
                    for beta = betas
                        obj.createFilterConnectivity(eta, beta, type);
                        [lambda, phi] = obj.computeEigenValueFunctional(n);
                        eigVs = [eigVs; lambda'];
                        eigFs = [eigFs; phi];
                    end
                    eigF1 = [eigF1; eigFs(:,1)'];
                    eigF2 = [eigF2; eigFs(:,2)'];
                    eigF3 = [eigF3; eigFs(:,3)'];
                    eigF4 = [eigF4; eigFs(:,4)'];
                    eigF5 = [eigF3; eigFs(:,5)'];
                    eigF6 = [eigF4; eigFs(:,6)'];
                    obj.plot(eigVs, betas, eta, n, analytical, type, j);
%                     obj.plotEig(eigFs, betas, eta, n, analytical, type, j);
                    j = j + 1;
                end
            end
            obj.createMonitoringEigenModes(1);
            obj.updateMonitoringEigenModes(eigF1);
            obj.createMonitoringEigenModes(2);
            obj.updateMonitoringEigenModes(eigF2);
            obj.createMonitoringEigenModes(3);
            obj.updateMonitoringEigenModes(eigF3);
            obj.createMonitoringEigenModes(4);
            obj.updateMonitoringEigenModes(eigF4);
            obj.createMonitoringEigenModes(5);
            obj.updateMonitoringEigenModes(eigF5);
            obj.createMonitoringEigenModes(6);
            obj.updateMonitoringEigenModes(eigF6);
        end

        function plot(obj, eigVs, betas, eta, n, analytical, type, j)
            subplot(2,3,j)
            for i=1:n
                plot(betas, eigVs(:,i),'.-')
                hold on;
            end
            for i=1:n
                plot(betas, ones(size(betas))*analytical(i),'--',Color='k')
                hold on;
            end
            legend([string(1:n),'Analytica  l'])
            title('\eta = '+string(eta)+', Type = '+string(type))
            xlabel('\beta') 
            ylabel('Eigenvalues') 
            ylim([-1 200])
            hold off;
            drawnow;
        end

        function updateMonitoringEigenModes(obj,eigFs) 
            for j = 1:size(eigFs,2)
                data = {};
                for i = 1:size(eigFs,1)
                    data = [data; [eigFs(i,j).fValues]];
                end
                obj.monitoringEigenModes.update(j,data);
                obj.monitoringEigenModes.refresh();
            end
        end

        function createMonitoringEigenModes(obj, j)
            nPlots = 6;
            chartTypes = cell(1,nPlots);
            barLims = cell(1,nPlots);
            funs = cell(1,nPlots);
            titles = cell(1,nPlots);
            field1 = obj.designVariable.fun;
            for i = 1:nPlots
                titles{i}      = string(i);
                chartTypes{i}   = 'surf';
                barLims{i} = [];
                funs{i} = field1;
            end
            s.shallDisplay = true;
            s.maxNColumns  = 3;
            etas = [0.0,0.5,1.0];
            s.titles       = ['Eig '+string(j)+', \eta = '+string(etas)+', Type = FP', 'Eig '+string(j)+' \eta = '+string(etas)+', Type = P'];
            s.chartTypes   = [chartTypes];
            s.barLims      = [barLims];
            s.funs         = [funs];
            obj.monitoringEigenModes = Monitoring(s);
        end 

        function createAnalyticalEigenModes(obj, m, n)
            x1 = linspace(0.25,0.75,100);
            x2 = linspace(0.25,0.75,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            mesh = Mesh.create(s);

            s.fHandle = @(x) cos(2*m*pi*(x(1,:,:)-0.25)).*cos(2*n*pi*(x(2,:,:)-0.25));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            fun       = AnalyticalFunction(s);
            fun.plot();
            drawnow;
        end           
    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.type        = 'RectangleInclusion';
            s.xSide       = 0.5;
            s.ySide       = 0.5;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
%             s.type        = 'ThreeRectanglesInclusion';
%             s.xSide1       = 0.3;
%             s.ySide1       = 0.3;
%             s.xCoorCenter1 = 0.4;
%             s.yCoorCenter1 = 0.5;
%             s.xSide2       = 0.05;
%             s.ySide2       = 0.05;
%             s.xCoorCenter2 = 0.6;
%             s.yCoorCenter2 = 0.5;
%             s.xSide3       = 1.0;
%             s.ySide3       = 0.1;
%             s.xCoorCenter3 = 0.5;
%             s.yCoorCenter3 = 0.2; 
            g             = GeometricalFunction(s);
            phi           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phi;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end        

        function createCharacteristicFunction(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);
            obj.characteristicFunction  = CharacteristicFunction.create(uMesh);
        end

        function createDesignVariable(obj)
            s.fun  = obj.filter.compute(obj.characteristicFunction,3);
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            dens.plot();
            obj.designVariable = dens;
        end

        function createFilterConnectivity(obj, eta, beta, type)
            if strcmp(type,"FP")
                s.filterType = 'FilterAndProject';
                s.mesh       = obj.mesh;
                s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%                 s.filterStep = 'LUMP';
                s.filterStep = 'PDE';
                s.beta       = beta;
                s.eta        = eta;
                f            = Filter.create(s);
                obj.filterConnect = f;
            elseif strcmp(type,"P")
                s.beta       = beta;
                s.eta        = eta;
                f = HeavisideProjector(s);
                obj.filterConnect = f;
            end
        end        

        function  [lambdas, phis] = computeEigenValueFunctional(obj, n)
            eigen = obj.computeEigenValueProblem();
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.mesh = obj.mesh;
            s.filter = obj.filterConnect;
            mE = MinimumEigenValueFunctional(s);
            [lambdas, phis] = mE.computeEigenModes(obj.designVariable, n);
        end

        function eigen = computeEigenValueProblem(obj)
            s.mesh  = obj.mesh;
            s.shift = 0.0;
            eigen   = StiffnessEigenModesComputer(s);
        end



    end

end
