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
            x = obj.xV;
            y = obj.yV;
            z = obj.value2print;
%             %figID = figure();
%             contour(x,y,z,50);
%             xlabel('$m_1$','Interpreter','latex');
%             ylabel('$m_2$','Interpreter','latex');
%             tN = obj.titleName;
%             title(['$',tN,'$'],'interpreter','latex');
%             colorbar;
%             
%             hold on                                   
%             v = [0,0];
%             [M,c] = contour(x,y,z,v);
%             c.LineWidth = 3;
%             c.LineColor = 'r';
                x2 = repmat(x,20,1);                    
                y2 = repmat(y',1,20);
%             
            s.fileName = fullfile(obj.outPutPath,[obj.fileName,'MxMy']);
            s.title    = obj.titleName;            
            s.axisAdder = MxMyAxisAdder();
            p =  SuperEllipseExponentContourPlotter(s);
            p.plot(x2(:),y2(:),z(:));            

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