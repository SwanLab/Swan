classdef StressNormPduringIterAsPostprocess < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        dataRes
        iterC
        cost
        fileNameIter
        finalIter
    end
    
    properties (Access = private)
        microCase
        path
        fileResName
        fileNameMesh
        pNorm
        valueT
        caseName
    end
    
    methods (Access = public)
            
        function obj = StressNormPduringIterAsPostprocess()
            obj.init();
            obj.obtainCost();
            value = zeros(length(obj.pNorm),obj.finalIter);
            for iter = 1:obj.finalIter
                fNameIter = [obj.fileResName,num2str(iter)];
                obj.wrapResAndMesh(fNameIter);
                value(:,iter) = obj.computeStressNorm();
                if mod(iter,1)== 0
                obj.plotValue(value,iter)
                end
            end
           obj.plotFinal(value);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.valueT = 0.18324;
            %obj.valueT  = 34.6295;
            %obj.valueT = 0.5646;
            %obj.microCase = 'Rectangle';
            obj.microCase = 'SuperEllipse';
            meshType = 'Small';
            %meshType  = 'Medium';
            obj.finalIter = 73;%607;%20;%46;%607;%42;18;405;42;296;57;46;250;485;250;
%             fCase = [obj.microCase,meshType];
%             obj.caseName = ['ExperimentingPlot'];
%             obj.path = ['/media/alex/My Passport/LatticeResults/StressNorm',fCase,'/'];
%             obj.fileResName  = ['ExperimentingPlot'];
            
            fCase = [obj.microCase,meshType];
            obj.caseName = ['ExperimentingPlot7',fCase];
            obj.path = ['/media/alex/My Passport/LatticeResults/',obj.caseName,'/'];
            obj.fileResName  = ['ExperimentingPlot',obj.microCase];
            obj.fileNameMesh = ['CantileverSquare',meshType];
            obj.pNorm = 2:2:16;
        end
        
        function obtainCost(obj)
            fNameMon = fullfile(obj.path,'Monitoring.fig');
            h = openfig(fNameMon,'invisible');
            handles = findobj(h,'Type','line');
            iter = get(handles(5),'Xdata');
            c = get(handles(5),'Ydata');
            close all
            obj.iterC = iter;
            obj.cost  = c;
        end
        
        function plotValue(obj,value,iter)
            figure(1)
            clf
            plot(obj.iterC,obj.cost*obj.valueT)
            hold on
            plot(1:iter,value(:,1:iter))   
            hold off
            drawnow
        end
        
        function plotFinal(obj,value)
            f = figure(1);
            clf
            p{1} = plot(obj.iterC,obj.cost*obj.valueT);
            hold on
            for i = 1:size(value,1)
              p{i+1} = plot(value(i,:));
            end
           pP = plotPrinter(f,p);
           fName = ['/home/alex/Dropbox/GregoireMeeting7Decembre/',obj.caseName];
           pP.print(fName)
        end
        
        function wrapResAndMesh(obj,fileNameIter)
            s.folderPath = obj.path;
            s.fileName   = fileNameIter;
            w = WrapperMshResFiles(s);
            w.compute();
            obj.mesh    = w.mesh;
            obj.dataRes = w.dataRes;
        end
        
         function value = computeStressNorm(obj)
            s.m1            = obj.dataRes.DesignVar1;
            s.m2            = obj.dataRes.DesignVar2;
            s.alpha         = obj.dataRes.AlphaGauss';
            switch obj.microCase 
                case 'Rectangle'
                    s.vademecumName = 'SuperEllipseQMax';
                case 'SuperEllipse'
                    s.vademecumName = 'SuperEllipseQOptAnalytic';
            end
            s.mesh          = obj.mesh;
            s.fileName      = obj.fileNameMesh;
            s.pNorm         = obj.pNorm;
            sC = StressNormFromVademecumComputer(s);
            value = sC.compute();
         end
    
    end
    
end