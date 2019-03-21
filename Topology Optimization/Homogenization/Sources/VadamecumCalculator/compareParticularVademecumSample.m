classdef compareParticularVademecumSample < handle
    
    properties (Access = private)
        pNorm
        fileName
        
        inclusionLength
        rectangleSmoothVM        
        rectangleVM
        monomials
    end
    
    methods (Access = public)
        
        function obj = compareParticularVademecumSample()
            obj.init();
            obj.computeSmoothRectangleSample();
            obj.computeRectangleSample();
            obj.plotAmplificators();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.pNorm = 8;
        end
        
        function computeSmoothRectangleSample(obj)
            obj.fileName = 'SmoothRectangle';
            obj.inclusionLength = 0.95;
            obj.rectangleSmoothVM = obj.computeCellVariables();            
        end
        
        function v = computeCellVariables(obj)
            obj.computeVademecum();
            v = obj.computePtensorVariables();
        end
        
        function computeVademecum(obj)
            d = obj.computeInputForVademecumCalculator();
            vc = VademecumCellVariablesCalculator(d);
            vc.computeVademecumData()
            vc.saveVademecumData();            
        end
        
        function vd = computePtensorVariables(obj)
            d.fileName = ['CaseOfStudy',obj.fileName];
            d.pNorm = obj.pNorm;
            vc = VademecumPtensorComputer(d);
            vc.compute();
            vd = vc.vademecumData;
        end
        
        function d = computeInputForVademecumCalculator(obj)
            d.fileName   = ['CaseOfStudy',obj.fileName];
            d.freeFemFileName = obj.fileName;
            d.mxMin = obj.inclusionLength;
            d.mxMax = obj.inclusionLength;
            d.myMin = obj.inclusionLength;
            d.myMax = obj.inclusionLength;
            d.nMx   = 1;
            d.nMy   = 1;
            d.outPutPath = [];            
        end
     
        function computeRectangleSample(obj)
            obj.fileName = 'Rectangle';
            obj.inclusionLength = obj.computeInclusionLengthForRectangle();
            obj.rectangleVM = obj.computeCellVariables();            
            obj.checkSameVolume();            
        end
        
        function m  = computeInclusionLengthForRectangle(obj)
            smoothVolume = obj.rectangleSmoothVM.variables{1,1}.volume;
            m = sqrt(1-smoothVolume);
        end
        
        function checkSameVolume(obj)
            sVolume = obj.rectangleSmoothVM.variables{1,1}.volume;
            rVolume = obj.rectangleVM.variables{1,1}.volume;
            isequal = (sVolume - rVolume)/rVolume < 1e-3;
            if ~isequal
                error('not same volume')
            end
        end
        
        function plotAmplificators(obj)
            y(:,1) = obj.rectangleSmoothVM.variables{1,1}.Ptensor;
            y(:,2) = obj.rectangleVM.variables{1,1}.Ptensor;
            fig = figure();
            h = bar(log(y));
            monomialsLeg = obj.computeAlphaLegend();
            set(gca, 'XTickLabel',monomialsLeg, 'XTick',1:numel(monomialsLeg));
            set(gca,'XTickLabelRotation',45);
            legend('SmoothRectangle','Rectangle')
            p = barPrinter(fig,h);
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            plotName = 'CaseOfStudy';
            p.print([path,plotName]);            
        end
        
       function a = computeAlphaLegend(obj)
            obj.monomials = obj.rectangleVM.monomials;
            for ia = 1:size(obj.monomials,1)
               a{ia} = mat2str(obj.monomials(ia,:)); 
            end
        end        
       
    end    
   
end