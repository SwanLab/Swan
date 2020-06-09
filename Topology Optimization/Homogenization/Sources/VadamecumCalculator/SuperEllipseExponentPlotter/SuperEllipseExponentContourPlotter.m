classdef SuperEllipseExponentContourPlotter < handle
    
    properties (Access = private)
        fig
    end
    
    properties (Access = private)
        fileName
        titleF
        axisAdder
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseExponentContourPlotter(cParams)
            obj.init(cParams);          
        end
        
        function plot(obj,x,y,z)
            obj.plotMesh(x,y,z)
            obj.addAxis();
            obj.addTitle();
            obj.print();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = [cParams.fileName,'Contour'];            
            obj.titleF   = cParams.title;
            obj.axisAdder = cParams.axisAdder;            
        end
        
        function plotMesh(obj,x,y,z)
            obj.fig = figure();                        
            connec = obj.obtainConnec(x,y);
            h = trisurf(connec,x,y,z);
            colorbar
            h.FaceColor = 'interp';
            h.EdgeColor = 'black';
            h.EdgeAlpha = 1;
            hold on
            plot3(x,y,10*z,'ro','MarkerFaceColor','red','MarkerSize',3)            
            view(2)              
        end

        function connec = obtainConnec(obj,x,y)
            connec = delaunay(x,y);  
            connec = obj.obtainQualityElements(connec,x,y);
        end
        
        function print(obj)
            fp = contourPrinter(obj.fig);
            filePath = obj.fileName;
            fp.print(filePath);             
        end
%         
%         function addAxis(obj)
%             xlabel(obj.xAxis,'Interpreter','latex');
%             set(gca,'xtick',[0:pi/8:pi/2]);
%             set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})                         
%             ylabel(obj.yAxis,'Interpreter','latex'); 
%         end
        
        function addAxis(obj)
            obj.axisAdder.add();
        end                  
        
        function addTitle(obj)
            tN = ['\textrm{',obj.titleF,'optimal super-ellipse exponent}'];
            title(['$',tN,'$'],'interpreter','latex')            
        end        
        
    end
    
    methods (Access = private, Static)
        
        

        function connec = obtainQualityElements(connec,x,y)            
            s.coord = [x,y];
            s.connec = connec;
            m = Mesh().create(s);
            qua = m.computeElementQuality';
            isQ = qua > 0.02;    
            connec = connec(isQ,:);
        end
        
    end
    
end