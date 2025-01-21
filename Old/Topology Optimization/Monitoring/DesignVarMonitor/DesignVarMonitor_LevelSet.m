classdef DesignVarMonitor_LevelSet < DesignVarMonitor_Abstract
    
    properties (Access = protected, Abstract)
        unfittedType
        meshIncludeBoxContour
    end
    
    properties (Access = protected)
        designVarName = 'Level Set - \phi';
        meshUnfitted
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet(cParams)
            obj@DesignVarMonitor_Abstract(cParams);
        end
        
        function plot(obj)
            obj.refreshFigure();
        end
        
    end
    
    methods (Access = protected)
        
        function initPlotting(obj)
             set(obj.axes,'CLim',[0, 1],'XTick',[],'YTick',[]);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [1 0 0];
        end
        
    end
    
    methods (Access = private)
        
        function refreshFigure(obj)
            figure(obj.figHandle.Number)
            cla reset;
            hold on
            uMesh = obj.designVar.getUnfittedMesh;
            uMesh.plot();
        end
        
    end
    
end