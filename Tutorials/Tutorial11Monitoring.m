classdef Tutorial11Monitoring < handle

    properties (Access = private)
        mesh1
        mesh2
        field1
        field2
        monitoring1
        monitoring2
    end 

    methods (Access = public)

        function obj = Tutorial11Monitoring()
            obj.init();
            % First example
            obj.createMonitoring1()
            obj.updateMonitoring1()
            % Second example
            obj.createMonitoring2()
            obj.updateMonitoring2()
        end

    end

    methods (Access = private)

        function init(obj)
            close all
            obj.createFields()
        end

        function createFields(obj)
            obj.mesh1 = QuadMesh(1,1,5,5);
            obj.mesh2 = TriangleMesh(1,1,5,5);
            obj.field1 = LagrangianFunction.create(obj.mesh1,1,'P1');
            obj.field2 = LagrangianFunction.create(obj.mesh2,2,'P1');
        end

        function createMonitoring1(obj)
            s.shallDisplay = true;
            s.maxNColumns = 2;
            s.titles = [{'Dataset X'},{'Dataset X&Y'},{'Datasets X'},{'Datasets X&Y'},{'SemiLogY'},{'LogLog'},{'Bar'}];
            s.chartTypes =[{'plot'},{'plot'},{'multiplot'},{'multiplot'},{'logy'},{'loglog'},{'bar'}]; 
            s.legends = [{["A","B","C"]},{["D","E","F"]}];
            obj.monitoring1 = Monitoring(s);
        end

        function updateMonitoring1(obj)
            for i=1:5
                PlotData1 = i;
                PlotData2 = [i;-i*5];
                MultiPlotData1 = [i;2*i;3*i];
                MultiPlotData2 = [i;2*i;3*i;10*i];
                BarData1 = [i;2*i];
                SemiLogData2 = [10*i;10*i];
                LogLogData = [i;10*i];
                obj.monitoring1.update(i,{[PlotData1],[PlotData2],[MultiPlotData1],[MultiPlotData2], ...
                                          [SemiLogData2],[LogLogData],[BarData1]})
                obj.monitoring1.refresh()
                pause(2);
            end
        end

        function createMonitoring2(obj)
            s.shallDisplay = true;
            s.maxNColumns = 3;
            s.titles = [{'Scalar Function'},{'Vector function X'},{'Vector function Y'}];
            s.chartTypes =[{'surf'},{'surf'},{'surf'}]; 
            s.funs = [{obj.field1},{obj.field2},{obj.field2}];
            s.barLims = [{[0 1]},{[]},{[]}];
            obj.monitoring2 = Monitoring(s);
        end

        function updateMonitoring2(obj)
            for i=1:5
                fv1 = (i/5)*ones(size(obj.field1.fValues));
                obj.field1.setFValues(fv1);

                fv2            = obj.field2.fValues;
                fv2(:,1)       = i;
                fv2(1:2:end,2) = 15-2.5*i;
                fv2(2:2:end,2) = -5+3.5*i;
                obj.field2.setFValues(fv2);

                SurfData1 = obj.field1.fValues;
                SurfData2 = obj.field2.fValues(:,1);
                SurfData3 = obj.field2.fValues(:,2);
                obj.monitoring2.update(i,{[SurfData1],[SurfData2],[SurfData3]})
                obj.monitoring2.refresh()
                pause(2);
            end 
        end

    end

end