classdef PlotterLevelSet < handle

    properties (Access = private)
        designVariable
        unfittedMesh
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
            %obj.unfittedMesh = cParams.unfittedMesh;
        end
        
        function createFigure(obj)
            obj.figHandle = figure();
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            obj.plotUnfittedMesh();
        end

        function plotUnfittedMesh(obj)
            hold on
            uMesh = obj.designVariable.getUnfittedMesh();
            %uMesh = obj.unfittedMesh;
            uMesh.plot();
        end
            
    end
    
end