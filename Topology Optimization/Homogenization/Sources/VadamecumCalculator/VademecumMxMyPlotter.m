classdef VademecumMxMyPlotter < VademecumPlotter
    
    properties (Access = protected)
        XYname = ' MxMy'
    end    
    
    methods (Access = public)
        
        function obj = VademecumMxMyPlotter(d)
            obj.init(d);
            obj.xV = obj.mxV;
            obj.yV = obj.myV;            
        end
        
        function plot(obj)
            obj.plotVolume();
            obj.plotHomogenizedTensor();
            obj.plotAmplificatorTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function plotFigure(obj)
            obj.fig = figure();
            x = obj.xV;
            y = obj.yV;
            z = obj.value2print;
            contour(x,y,z,50);
            xlabel('mx');
            ylabel('my');
            obj.addTitle();
            colorbar;
        end
        
    end
    
    methods (Access = private)
        
        function plotVolume(obj)
            obj.fileName = 'Volume';
            obj.value2print = obj.volume;
            obj.plotFigure();
            obj.printFigure();
        end       
        
    end    
    
end