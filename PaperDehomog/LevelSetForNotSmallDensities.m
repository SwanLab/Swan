classdef LevelSetForNotSmallDensities < handle

    properties (Access = private)      
        experimentData
        mesh
        density
    end

    methods (Access = public)

        function obj = LevelSetForNotSmallDensities()
            obj.loadDataExperiment();
            obj.createMesh();
            obj.createDensity();
            obj.createLevelSetFunction();
        end

    end

    methods (Access = private)

        function loadDataExperiment(obj)
           d = load('DataExampleArchCoarse.mat');  
           w = d.w;
           obj.experimentData = w;
        end

        function createMesh(obj)
            d = obj.experimentData;
            s.coord  = d.mesh.coord;
            s.connec = d.mesh.connec;
            obj.mesh = Mesh.create(s);
        end     

        function createDensity(obj)
            rhoV = obj.experimentData.dataRes.DensityGauss;
            rho  = LagrangianFunction.create(obj.mesh,1,'P0');
            rho.setFValues(rhoV);
            obj.density = rho;
            obj.density.plot
        end        

        function createLevelSetFunction(obj)

            %x1 = @(x) x(1,:,:);
            %x2 = @(x) x(2,:,:);
            
           % l  = 1;
           % x0 = 0;
           % y0 = 0;
           % fH = @(x) max(abs(x1(x)-x0),abs(x2(x)-y0))/l - 0.5;
            fH = 1-2*Heaviside(obj.density - 0.05);
            fH = project(fH,'P1');

            levelSet  = fH;
            
            figure()
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet.fValues);            
            uMesh.plot

        end

    end

end