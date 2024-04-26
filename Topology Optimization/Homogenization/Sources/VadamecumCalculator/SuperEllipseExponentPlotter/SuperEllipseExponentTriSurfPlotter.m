classdef SuperEllipseExponentTriSurfPlotter < handle
    
    properties  (Access = private)
        fig
    end
    
    properties (Access = private)
        fileName
        axisAdder
        titleF
        addAxisFunc
        viewV
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseExponentTriSurfPlotter(cParams)
            obj.init(cParams);          
        end
        
        function plot(obj,x,y,z)
            obj.plotMesh(x,y,z);
            obj.addTitle();            
            obj.addAxis();
            obj.print();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = [cParams.fileName,'TriSurf'];
            obj.titleF = cParams.title;    
            obj.axisAdder = cParams.axisAdder;      
            obj.viewV = cParams.view;
        end
        
        function plotMesh(obj,x,y,z)           
            obj.fig = figure();
            connec = obj.obtainConnec(x,y);                        
            trisurf(connec,x,y,z);
            hold on
            plot3(x,y,z,'ro','MarkerFaceColor','red','MarkerSize',3)
            view(obj.viewV(1),obj.viewV(2))
        end
        
        function connec = obtainConnec(obj,x,y)
            connec = delaunay(x,y);  
            connec = obj.obtainQualityElements(connec,x,y);
        end        
        
        function print(obj)
            fp = surfPrinter(obj.fig);
            filePath = obj.fileName;
            fp.print(filePath);             
        end
        
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
            m = Mesh.create(s);
            qua = m.computeElementQuality';
            isQ = qua > 0.02;    
            connec = connec(isQ,:);
        end
        
    end    
       
end