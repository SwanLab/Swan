classdef OptimalSuperEllipseExponentPlotter < handle
    
    properties (Access = private)
        outputPath
        averagedSuperEllipse
        vademecum
        fig
    end
    
    methods (Access = public)
        
        function obj = OptimalSuperEllipseExponentPlotter()
            obj.outputPath = '/home/alex/git-repos/MicroStructurePaper/';
            obj.readVademecum();
            obj.computePonderatedOptimalSuperEllipse();
            obj.plotAverage();
            obj.plotDeviation();
          end
        
    end
    
    methods (Access = private)
        
        function readVademecum(obj)
            v = VademecumReader();   
            obj.vademecum = v;         
        end
        
        function computePonderatedOptimalSuperEllipse(obj)
            s.vademecum = obj.vademecum;
            p = PonderatedOptimalSuperEllipseComputer(s);
            p.compute();
            obj.averagedSuperEllipse = p;
        end
        
        function plotAverage(obj)
            s.xiV = obj.vademecum.xiV;
            s.rhoV = obj.vademecum.rhoV;
            s.qMean = obj.averagedSuperEllipse.qMean;                
            s.value = obj.averagedSuperEllipse.qMean;    
            s.title = 'Average of ';
            s.fileName = [obj.outputPath,'Qmean'];
            p = SuperEllipseExponentPlotter(s);            
            p.plot();            
        end
        
        function plotDeviation(obj)
            s.xiV = obj.vademecum.xiV;
            s.rhoV = obj.vademecum.rhoV;
            s.qMean = obj.averagedSuperEllipse.qMean;            
            s.value = obj.averagedSuperEllipse.qDesv;  
            s.title = 'Deviation of ';            
            s.fileName = [obj.outputPath,'Qdesv'];
            p = SuperEllipseExponentPlotter(s);            
            p.plot();            
        end
                
    end    

end