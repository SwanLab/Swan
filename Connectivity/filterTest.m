classdef filterTest < handle
    
    properties (Access = public)
        mesh
        levelSet
        characteristicFunction
        filter
        filterAndProject
        projector
        density
        densityNR
        densityR
        densityF
        densityFP
        densityP
    end

    methods (Access = public)

        function obj = filterTest()
            close all;
            obj.createMesh()
            obj.createDensity()
            obj.createFilterAndProject()
            obj.createFilter()
            obj.createHeavisideProjector()
            obj.createFilterAndProject()
            obj.computeDensity()
        end
        
        function createMesh(obj)
            x1 = linspace(0,2,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createDensity(obj)
            s.fHandle = @(x) x(1,:,:)/2;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            obj.density        = DesignVariable.create(sD);
            obj.density.plot()
            obj.density.fun.plot()

        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createHeavisideProjector(obj)
            s.beta = 10.0;
            s.eta = 0.5;
            obj.projector = HeavisideProjector(s);
        end 

        function createFilterAndProject(obj)
            s.beta = 16.0;
            s.eta = 0.0;
            s.filterStep = 'LUMP';
            s.mesh = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.filterAndProject = FilterAndProject(s);
        end 


        function computeDensity(obj, x)          
            % Not Rounding Densities
            s.operation = @(xV) obj.computeComplementaryDensity(obj.density.fun,xV);
            obj.densityNR = DomainFunction(s);
            obj.densityNR.project('P1',obj.mesh).plot()

            % Rounding Densities
            s.operation = @(xV) obj.computeRoundedComplementaryDensity(obj.density.fun,xV);
            obj.densityR = DomainFunction(s);
            obj.densityR.project('P1',obj.mesh).plot()

            % Filter
            obj.densityF = obj.filter.compute(1-obj.density.fun, 2);
            obj.densityF.plot()

            % Project
            fV = obj.projector.project(1-obj.density.fun);
            obj.densityP = LagrangianFunction.create(obj.density.fun.mesh,1,'P1')
            obj.densityP.fValues = fV
            obj.densityP.plot()

            % Filter and Project
            obj.densityFP = obj.filterAndProject.compute(1-obj.density.fun, 2);
            obj.densityFP.plot()
        end
        
        function rho = computeComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
        end

        function rho = computeRoundedComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
            rho = round(rho);
            rho = max(0,min(1,rho));
         end


    end
end