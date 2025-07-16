classdef PlotterLevelSet < handle

    properties (Access = private)
        designVariable
    end

    properties (Access = private)
        figHandle
    end

    methods (Access = public)

        function obj = PlotterLevelSet(cParams)
            obj.init(cParams);
            obj.createFigure();
        end
        
        function plot(obj)
            figure(obj.figHandle);
            cla reset;
            obj.plotUnfittedMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end
        
        function createFigure(obj)
            obj.figHandle = figure('units','normalized','outerposition',[0 0 1 1]); %figure();
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            obj.plotUnfittedMesh();
        end

        function plotUnfittedMesh(obj)
            hold on
            uMesh = obj.designVariable.getUnfittedMesh();
            uMesh.plot();
        end
            
    end
    
end