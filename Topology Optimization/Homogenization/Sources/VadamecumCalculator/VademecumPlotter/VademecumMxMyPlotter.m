classdef VademecumMxMyPlotter < VademecumPlotter
    
    properties (Access = protected)
       name = 'MxMy'; 
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
            obj.plotHomogenizedTensorIsotropy();
            obj.plotAmplificatorTensor();
        end
        
    end
    
    methods (Access = protected)
        
        function plotFigure(obj)
            obj.fig = figure();
            x = obj.xV;
            y = obj.yV;
            z = obj.value2print;
            %figID = figure();
            contour(x,y,z,50);
            xlabel('$m_1$','Interpreter','latex');
            ylabel('$m_2$','Interpreter','latex');
            tN = obj.titleName;
            title(['$',tN,'$'],'interpreter','latex');
            colorbar;
            
            hold on                                   
            v = [0,0];
            [M,c] = contour(x,y,z,v);
            c.LineWidth = 3;
            c.LineColor = 'r';
            
%             s.x = x;
%             s.y = y;
%             s.z = z;  
%             s.figID = figID;
%             lineAdder = ZeroLinePlotAdder(s);
%             lineAdder.addLine();            
        end
        
    end
    
    methods (Access = private)
        
        function plotVolume(obj)
            obj.fileName = 'Volume';
            obj.titleName = '\textrm{Volume}';
            obj.value2print = obj.volume;
            obj.plotFigure();
            obj.printFigure();
        end       
        
    end    
    
   
    
end