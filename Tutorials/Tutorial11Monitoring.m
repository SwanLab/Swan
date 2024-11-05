classdef Tutorial11Monitoring < handle

    properties (Access = private)
        mesh
        field
        monitoring
    end

    methods (Access = public)

        function obj = Tutorial11Monitoring()
            obj.init();
            obj.createMonitoring()
            obj.updateMonitoring()
        end

    end

    methods (Access = private)

        function init(obj)
            close all
            obj.createField()
        end

        function createField(obj)
            obj.mesh = QuadMesh(1,1,5,5);
            obj.field = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.field.fValues(:) = 1;
        end

        function createMonitoring(obj)
            s.shallDisplay = true;
            s.maxNColumns = 3;
            s.titles = [{'counting'},{'counting multiples of 5'},{'field'}];%,{'counting & doubles'},{'counting & doubles of multiples of 5'}
            s.chartTypes =[{'plot'},{'surf'},{'plot'}]; %,{'multiPlot'},{'multiPlot'}
            s.mesh = obj.mesh;
            s.barLim = [0 10];
            s.legend = ["A","B","C"];
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            for i=1:5
                PlotData1 = [1*i 2*i 3*i i];
                PlotData2 = [1 2 3];
                SurfData(:) = 2*i.*ones(size(obj.field.fValues));
                obj.monitoring.update(i,{[PlotData1],[SurfData],[PlotData2]})
                obj.monitoring.refresh()
                pause(2);
            end
        end

    end

end