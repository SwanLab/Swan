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
            s.titles = [{'counting'},{'field'},{'counting multiples of 5'}];%,{'counting & doubles'},{'counting & doubles of multiples of 5'}
            s.chartTypes =[{'plot'},{'surf'},{'plot'}]; %,{'multiPlot'},{'multiPlot'}
            s.mesh = obj.mesh;
            s.barLim = [0 10];
            s.legend = ["A","B","C"];
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            for i=1:5
                PlotData1 = i;
                PlotData2 = [i;i*5];
                MultiPlotData = [i;i*5];
                MultiPlotData2 = [i;i*5;i*10];
                SurfData(:) = 2*i.*ones(size(obj.field.fValues));
                SurfData(1:6) = 3*i;
                obj.monitoring.update(i,{[PlotData1],[SurfData],[PlotData2]})
                obj.monitoring.refresh()
                pause(2);
            end
        end

    end

end