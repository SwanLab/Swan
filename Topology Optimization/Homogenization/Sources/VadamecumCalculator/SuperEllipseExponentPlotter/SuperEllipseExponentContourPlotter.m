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
        
        function print(obj)
            fp = contourPrinter(obj.fig);
            filePath = obj.fileName;
            fp.print(filePath);             
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
            h.EdgeAlpha = 0.3;
            h.FaceAlpha = 1;            
            hold on
            plot3(x,y,10*abs(z),'ro','MarkerFaceColor','red','MarkerSize',3)            
            view(2)              
        end

        function connec = obtainConnec(obj,x,y)
            connec = delaunay(x,y);  
            connec = obj.obtainQualityElements(connec,x,y);
        end
        
        function addAxis(obj)
            obj.axisAdder.add();
        end                  
        
        function addTitle(obj)
            tN = ['\textrm{',obj.titleF,'optimal super-ellipse exponent}'];
          %  title(['$',tN,'$'],'interpreter','latex')            
        end        
        
    end
    
    methods (Access = private, Static)
     
        function newConnec = obtainQualityElements(connec,x,y)            
            s.x = x;
            s.y = y;
            s.connec = connec;
            c = ConnecWithQualityComputer(s);
            newConnec = c.compute();
        end
        
    end
    
end