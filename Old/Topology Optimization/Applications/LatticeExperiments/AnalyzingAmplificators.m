classdef AnalyzingAmplificators < handle
    
    properties (Access = private)
       homog
       m1
       m2
    end
    
    properties (Access = private)
        fileName
        pNorm
        path
    end
    
    
    methods (Access = public)
        
        function obj = AnalyzingAmplificators()
            obj.init();
            obj.createM1M2();
            obj.createHomogenizedComputer();
            obj.computeAmplificators()
            obj.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'SuperEllipseQOptAnalytic';
            obj.pNorm    = 16;
            obj.path     = '/home/alex/Dropbox/GregoireMeeting19January/';
        end
        
        function createM1M2(obj)
            npoints = 100;
            obj.m1 = 0.5*ones(npoints,1);
            obj.m2 = linspace(0.8,0.97,npoints)';
        end
        
        function createHomogenizedComputer(obj)
           s.type = 'ByVademecum';
           s.fileName = obj.fileName;
           h = HomogenizedVarComputer.create(s);
           obj.homog = h;
        end
        
        function computeAmplificators(obj)
            x{1} = obj.m1;
            x{2} = obj.m2;
            x{3} = [ones(1,length(obj.m1));zeros(1,length(obj.m1))];
            obj.homog.computePtensor(x,obj.pNorm);
        end
        
        function plot(obj)
            obj.plot1111();
            obj.plot2222();
            obj.plot1212();
            obj.printFigure(f,h)
        end
        
        function plot1212(obj)
            comp = 165;
            name = 'P1212';
            obj.plotComponent(comp,name);
        end
        
        function plot1111(obj)
            comp = 1287;
            name = 'P1111';
            obj.plotComponent(comp,name);
        end
        
        function plot2222(obj)
            comp = 495;
            name = 'P2222';
            obj.plotComponent(comp,name);
        end
        
        function plotComponent(obj,comp,name)
            f = figure();
            Pp = obj.homog.Pp(comp,:);
            Ppn = (Pp).^(1/obj.pNorm);
            h{1} = semilogy(obj.m2,(Ppn),'+-');
            legend(name);
            obj.printFigure(f,h,name)
        end
        
        function printFigure(obj,f,h,name)
            outputName = [obj.path,name,'m2_05'];
            printer = plotPrinter(f,h);
            printer.print(outputName);
        end
        
    end
    
end