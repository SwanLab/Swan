classdef ConnectivityComputer < handle

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

        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable();  
            obj.createFilterConnectivity();
            [eigV, eigF] = obj.computeEigenValueFunctional();
            obj.createMonitoringEigenModes(eigV);
            obj.updateMonitoringEigenModes(eigF);
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
%             s.fValues = round(obj.characteristicFunction.project('P1').fValues);
%             s.mesh    = obj.mesh;
%             s.order   = 'P1';
%             s.fun = LagrangianFunction(s);
            s.fun  = obj.filter.compute(obj.characteristicFunction,3);
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            dens.plot();
            obj.designVariable = dens;
        end

        function createFilterConnectivity(obj)
%            s.filterType = 'FilterAndProject';
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'LUMP';
%             s.beta       = 100.0;
%             s.eta        = 1.0;
%             f            = Filter.create(s);
%             obj.filterConnect = f;
%             s.filterType = 'FilterAdjointAndProject';
%             f            = Filter.create(s);
%             obj.filterAdjointConnect = f;
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
%             f.updateEpsilon(1.0*obj   .mesh.computeMeanCellSize());
            obj.filterConnect = f;

        end        

        function [eigV, eigF] = computeEigenValueFunctional(obj)
            eigen = obj.computeEigenValueProblem();
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.mesh = obj.mesh;
            s.filter = obj.filterConnect;
%             s.filterAdjoint = obj.filterAdjointConnect;
            mE = MinimumEigenValueFunctional(s);
%             [lambda, dlambda] = mE.computeFunctionAndGradient({obj.designVariable,0});  
            [eigV, eigF] = mE.computeEigenModes(obj.designVariable, 6);
        end

        function eigen = computeEigenValueProblem(obj)
            s.mesh  = obj.mesh;
            s.shift = 0.0;
            eigen   = StiffnessEigenModesComputer(s);
        end

        function createMonitoringEigenModes(obj,eigV)
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
            s.titles       = ['Eig '+string(1:6)+'= '+string(eigV)'];
            s.chartTypes   = [chartTypes];
            s.barLims      = [barLims];
            s.funs         = [funs];
            obj.monitoringEigenModes = Monitoring(s);
        end 


        function updateMonitoringEigenModes(obj,eigF) 
            data = {};
            for i = 1:size(eigF,2)
                data = [data; [eigF(i).fValues]];
            end
            obj.monitoringEigenModes.update(1,data);
            obj.monitoringEigenModes.refresh();
        end

    end

end