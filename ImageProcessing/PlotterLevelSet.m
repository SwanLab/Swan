classdef PlotterLevelSet < handle

    properties (Access = private)
        designVariable
    end

    properties (Access = private)
        figHandle
    end

    methods (Access = public)

        function obj = PlotterLevelSet(cParams)
            obj.createFigure(cParams.unfittedMesh);
        end
        
        function plot(obj,uMesh)
            figure(obj.figHandle);
            cla reset;
            hold on
            uMesh.plot();
        end
        
    end
    
    methods (Access = private)
        
        
        function createFigure(obj,uMesh)
            obj.figHandle = figure();
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            uMesh.plot();            
        end


    end
    
end