classdef testAmplificatorTensorForInclusions < ...
        testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = private)
        ampTensorNum        
        fileName
        printTopology  
        m2
        m1
        PinR
        CR
        volR
        PinS
        CS
        volS
    end
    
    properties (Access = protected)
        variablesToStore = {'P'};
        tol = 1e-6;
        testName = 'AmplificatorTensorForInclusion';
    end
    
    methods (Access = public)
        
        function obj = testAmplificatorTensorForInclusions()
            obj.createNumericalAmplificationTensor()
            obj.selectComputedVar()            
        end
        
    end
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.computedVar{1} = obj.ampTensorNum.getValue();
        end
        
    end
    
    methods (Access = private)
        
        function createNumericalAmplificationTensor(obj)
            obj.m2 = 0.2;
            obj.printTopology  = true;
            obj.fileName       = strcat(obj.testName);
            n = 15;
            obj.PinR = zeros(9,n);
            obj.CR   = zeros(9,n);
            obj.volR = zeros(n,1);
            obj.PinS = zeros(9,n);
            obj.CS   = zeros(9,n);
            obj.volS = zeros(n,1);

            m1vect = linspace(0.001,0.99,n);
            for im1 = 1:n
                obj.m1 = m1vect(im1);
                obj.computeVariablesForRectangularInclusion(im1);
                obj.computeVariablesForSmoothRectangularInclusion(im1);
            end
            obj.plotComparison()
            obj.ampTensorNum = SymmetricFourthOrderPlaneStressVoigtTensor();
            obj.ampTensorNum.createRandomTensor();
        end   
        
        
        function plotComparison(obj)
            indexToPlot = [1 2 3 5 6 9];
            legP = {'P11', 'P12','P13','P22','P23','P33'};
            legC = {'C11', 'C12','C13','C22','C23','C33'};
            for i = 1:length(indexToPlot)
                figureID = figure(i);
                in = indexToPlot(i);
                h1 = plot(obj.volR,[obj.CR(in,:)',obj.PinR(in,:)']','-+');
                hold on
                h2 = plot(obj.volS,[obj.CS(in,:)',obj.PinS(in,:)']','-+');
                lege{1} = [legC{i},' rect inclusion '];
                lege{2} = [legP{i},' rect inclusion '];
                lege{3} = [legC{i},' smooth rect inclusion '];
                lege{4} = [legP{i},' smooth rect inclusion '];
                legend(lege{:},'Location','Best');
                xlabel('volume');
                obj.printFigure(figureID,h1,h2,legP{i});
            end
        end
        
        
        function printFigure(obj,figureID,h1,h2,lege)            
            outPutFile = ['/home/alex/Dropbox/Amplificators/Images/',lege];
            h{1} = h1;
            h{2} = h2;
            fp = plotPrinter(figureID,h);
            fp.print(outPutFile);
        end
        
        function computeVariablesForRectangularInclusion(obj,im1)
            homog              = obj.createHomogenizerWithRectangularInclusion(im1);
            Ch                 = homog.getCh();
            Ptensor            = homog.getAmplificatorTensor();
            PTensorValues      = Ptensor.getValue();
            Pinv               = obj.computePinvTensor(PTensorValues);
            obj.PinR(:,im1)    = Pinv(:);
            obj.CR(:,im1)      = Ch(:);
            obj.volR(im1)      = homog.getVolume();            
        end
        
        function computeVariablesForSmoothRectangularInclusion(obj,im1)
            homog              = obj.createHomogenizerWithSmoothRectangularInclusion(im1);
            Ch                 = homog.getCh();
            Ptensor            = homog.getAmplificatorTensor();
            PTensorValues      = Ptensor.getValue();
            Pinv               = obj.computePinvTensor(PTensorValues);
            obj.PinS(:,im1)    = Pinv(:);
            obj.CS(:,im1)      = Ch(:);
            obj.volS(im1)      = homog.getVolume();
        end
        
        function Pinv = computePinvTensor(obj,P)
            a = CompliancePlaneStressVoigtTensor;
            a.setValue(P);
            b = Inverter.invert(a);
            Pinv = b.getValue;
        end
        
        function h = createHomogenizerWithRectangularInclusion(obj,iter)
            f     = obj.fileName;
            m1v   = obj.m1;
            m2v   = obj.m2;
            print = obj.printTopology;
            h     = NumericalRectangleHomogenizer(f,print,m1v,m2v,iter);
        end
        
        function h = createHomogenizerWithSmoothRectangularInclusion(obj,iter)
            f     = obj.fileName;
            m1v   = obj.m1;
            m2v   = obj.m2; 
            print = obj.printTopology;
            h     = NumericalSmoothRectangleHomogenizer(f,print,m1v,m2v,iter);
        end
        
    end

end

