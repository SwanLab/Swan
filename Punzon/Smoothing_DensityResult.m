classdef Smoothing_DensityResult < handle

    properties (Access = private)
        mesh
        filter
        designVariable
    end

    methods (Access = public)

        function obj = Smoothing_DensityResult()
            obj.init()
            obj.createMesh();
            obj.createFilter();
            obj.createDesignVariable();
        end

    end

    methods (Access = private)



        function createMesh(obj)
            file = 'punzon3';
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            % ls= obj.computeLevelSet();
            % sUm.backgroundMesh = obj.mesh;
            % sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            % uMesh              = UnfittedMesh(sUm);
            % uMesh.compute(ls);
            % 
            % funLS        = CharacteristicFunction.create(uMesh);
            % f = obj.filter;
            % s.fun      = f.compute(funLS,2);
            % s.type     = 'Density';
            % s.plotting = true;
            % dens        = DesignVariable.create(s);
            % obj.designVariable = dens;
            load("densityProva.mat");
            s.fValues=density;
            s.mesh=obj.mesh;
            s.order='P1';
            f=LagrangianFunction(s);
            ls2 = f-0.5;
            ls2=ls2.project('P1');
            
            for i=1:30
                ls2.print(["provaa"+num2str(i)]);
                ls2=obj.filter.compute(ls2,2);
            end

        end


        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function ls=computeLevelSet(obj)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0.5;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = 0.4;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(obj.mesh);
            ls        = lsFun.fValues;
        end


    end

    methods(Static,Access=private)
        function init()
            close all;
        end  
    end
end
