classdef compareAmplificatorsWithMeshSize < handle
    
    properties (Access = private)
        pNorm
        fullFileName
        fileName
        hMax
        hMaxValues
        inclusionLength
        
        iter
        rectangleVM
        monomials
        Pcomp
        legends
    end
    
    methods (Access = public)
        
        function obj = compareAmplificatorsWithMeshSize()
            obj.init();
            obj.computeSamplesForDifferentMeshes();
            obj.plotAmplificators();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.pNorm = 8;
            obj.fileName = 'Rectangle';
            obj.inclusionLength = 0.913525506834508;%0.95;
            obj.hMaxValues = [0.04 0.02 0.01 0.005 0.0025];
        end
        
        function computeSamplesForDifferentMeshes(obj)
            for ih = 1:length(obj.hMaxValues)
                obj.iter = ih;
                obj.hMax = obj.hMaxValues(ih);
                obj.computeSampleForParticularMesh();
            end
        end
        
        function computeSampleForParticularMesh(obj)
            obj.fullFileName = ['CaseOfStudy',obj.fileName,num2str(obj.iter)];
            obj.rectangleVM{obj.iter} = obj.computeCellVariables();
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
            d.fileName = obj.fullFileName;
            d.pNorm = obj.pNorm;
            vc = VademecumPtensorComputer(d);
            vc.compute();
            vd = vc.vademecumData;
        end
        
        function d = computeInputForVademecumCalculator(obj)
            d = SettingsVademecumCellVariablesCalculator();
            d.fileName   = obj.fullFileName;
            d.freeFemFileName = obj.fileName;
            d.mxMin = obj.inclusionLength;
            d.mxMax = obj.inclusionLength;
            d.myMin = obj.inclusionLength;
            d.myMax = obj.inclusionLength;
            d.nMx   = 1;
            d.nMy   = 1;
            d.freeFemSettings.hMax = obj.hMax;
            d.outPutPath = [];
        end
        
        function plotAmplificators(obj)
            obj.obtainPcomponents();
            fig = figure();
            h = bar(log10(obj.Pcomp));
            monomialsLeg = obj.computeAlphaLegend();
            set(gca, 'XTickLabel',monomialsLeg, 'XTick',1:numel(monomialsLeg));
            set(gca,'XTickLabelRotation',45);
            obj.createLegend();
            legend(obj.legends,'Interpreter','latex')
            p = barPrinter(fig,h);
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            plotName = 'CaseOfStudyMeshVariation';
            p.print([path,plotName]);
        end
        
        function createLegend(obj)
            for ih = 1:length(obj.hMaxValues)
                obj.legends{ih} = ['$h=',num2str(obj.hMaxValues(ih)),'$'];
            end
        end
        
        function obtainPcomponents(obj)
            for ih = 1:length(obj.hMaxValues)
                comp = obj.rectangleVM{ih}.variables{1,1}.Ptensor;
                obj.Pcomp(:,ih) = comp;
            end
        end
        
        function a = computeAlphaLegend(obj)
            obj.monomials = obj.rectangleVM{1}.monomials;
            for ia = 1:size(obj.monomials,1)
                a{ia} = mat2str(obj.monomials(ia,:));
            end
        end
        
    end
    
end