classdef compareParticularVademecumSample < handle
    
    properties (Access = private)
        pNorms
        pNorm
        qNorms
        qNorm
        fileName
        volume
        volumes
        inclutionRatio
        
        
        qStrValue

        vademecums
        monomials
        prefixName
        nonShearIndex
        Pcomp
        print
        legends
        plottingShearTerms
        incLength
    end
    
    methods (Access = public)
        
        function obj = compareParticularVademecumSample()
            obj.init();
            obj.computeVademecum();
            obj.computeAndPlotPtensors();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.volumes = 0.95;%0.95;%[0.7,0.95];%0.1;
            obj.inclutionRatio = 1;%1;%0.5;
            obj.pNorms = [2,4,8];
            obj.qNorms = [2,4,8,16,32];%[2,4,8,16,32];%[4,8,16,32];
            obj.plottingShearTerms = true;
        end
        
        function computeAndPlotPtensors(obj)
            for ivolume = 1:length(obj.volumes)
                obj.volume = obj.volumes(ivolume);
                for ip = 1:length(obj.pNorms)
                    obj.pNorm = obj.pNorms(ip);
                    for iq = 1:length(obj.qNorms)
                        obj.qNorm = obj.qNorms(iq);
                        obj.computeSmoothRectanglePtensorSamples(iq);
                    end
                    obj.computeRectanglePtensorSamples();
                    obj.plotAmplificators();
                end
            end
        end
        
        function computeVademecum(obj)
            for ivolume = 1:length(obj.volumes)
                obj.volume = obj.volumes(ivolume);
                obj.computeSmoothRectangleSamples()
                obj.computeRectangleSample();
            end
        end
        
        function computeSmoothRectangleSamples(obj)
            for iq = 1:length(obj.qNorms)
                obj.qNorm = obj.qNorms(iq);
                obj.computeSmoothRectangleSample();
            end
        end
        
        function computeSmoothRectangleSample(obj)
            d.vademecumCase  = 'SmoothRectangle';
            obj.computeSample(d)
        end
        
        function computeRectangleSample(obj)
            d.vademecumCase  = 'Rectangle';
            obj.computeSample(d)
        end
        
        function computeSample(obj,d)
            d.volume         = obj.volume;
            d.inclutionRatio = obj.inclutionRatio;
            d.qNorm          = obj.qNorm;
            v = VademecumComputerForGivenVolume.create(d);
            v.compute();
        end        
        
        function computeSmoothRectanglePtensorSamples(obj,iq)
            obj.fileName = 'SmoothRectangle';
            v = obj.computeSmoothPtensorSample();
            obj.vademecums{iq,1} = v;
        end
        
        function computeRectanglePtensorSamples(obj)
            obj.fileName = 'Rectangle';
            v = obj.computeRectanglePtensorSample();
            nq = length(obj.qNorms);
            obj.vademecums{nq+1,1} = v;
        end
        
        function v = computeSmoothPtensorSample(obj)
            obj.fileName  = 'SmoothRectangle';
            obj.qStrValue = obj.qNorm;
            v = obj.computePtensorSample();
        end
        
        function v = computeRectanglePtensorSample(obj)
            obj.fileName  = 'Rectangle';
            obj.qStrValue = Inf;
            v = obj.computePtensorSample();
        end
        
        function v = computePtensorSample(obj)
            obj.obtainPrefixName();
            obj.print = true;
            v = obj.computePtensorVariables();
        end
        
        function obtainPrefixName(obj)
            volumeStr = obj.num2strIncludingDots(obj.volume);
            txiStr    = obj.num2strIncludingDots(obj.inclutionRatio);
            qStr      = obj.num2strIncludingDots(obj.qStrValue);                          
            obj.prefixName = ['CaseOfStudy','Rho',volumeStr,'Txi',txiStr,'Q',qStr];
        end
        
        function vd = computePtensorVariables(obj)
            d.fileName = [obj.prefixName,obj.fileName];
            d.pNorm = obj.pNorm;
            vc = VademecumPtensorComputer(d);
            vc.compute();
            vd = vc.vademecumData;
        end        
        
        function plotAmplificators(obj)
            obj.obtainMonomials();
            obj.obtainNonShearTerms();
            obj.obtainPcomponents();
            obj.makePlot();
        end
        
        function makePlot(obj)
            fig = figure();
            y = obj.Pcomp(obj.nonShearIndex,:);
            ind = y<0;
            y = (abs(y).^(1/obj.pNorm));
            y(ind) = -y(ind);
            h = bar(y);
            allMonLegends = obj.computeAlphaLegend();
            monomialsLeg = allMonLegends(obj.nonShearIndex);
            set(gca, 'XTickLabel',monomialsLeg, 'XTick',1:numel(monomialsLeg));
            set(gca,'XTickLabelRotation',45);
            obj.createLegend();
            legend(obj.legends,'Interpreter','latex','Location','Best')
            p = barPrinter(fig,h);
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            volumeStr = strrep(num2str(obj.volume),'.','');
            plotName = ['Rho',volumeStr,'P',num2str(obj.pNorm)];
            if obj.plottingShearTerms
                plotName = [plotName,'WithShear'];
            end
            p.print([path,plotName]);
        end
        
        function obtainPcomponents(obj)
            nvad = size(obj.vademecums,1);
            obj.Pcomp = zeros(size(obj.monomials,1),nvad);
            for ivad = 1:nvad
                Ptensor = obj.vademecums{ivad,1}.variables{1,1}.Ptensor;
                obj.Pcomp(:,ivad) = Ptensor;
            end
        end
        
        function a = computeAlphaLegend(obj)
            for ia = 1:size(obj.monomials,1)
                a{ia} = mat2str(obj.monomials(ia,:));
            end
        end
        
        function obtainMonomials(obj)
            obj.monomials = obj.vademecums{1}.monomials;
        end
        
        function obtainNonShearTerms(obj)
            nMon = size(obj.monomials,1);
            obj.nonShearIndex = true(nMon,1);
            if ~obj.plottingShearTerms
                for ia = 1:nMon
                    mon = obj.monomials(ia,:);
                    obj.nonShearIndex(ia) = obj.hasNoShearComponent(mon);
                end
            end
        end
        
        function createLegend(obj)
            nQnorms = length(obj.qNorms);
            for ih = 1:nQnorms
                obj.legends{ih} = ['$q=',num2str(obj.qNorms(ih)),'$'];
            end
            obj.legends{nQnorms+1} = 'Rectangle';
        end
        
    end
    
    methods (Access = private, Static)
        
        function itHasNoShear = hasNoShearComponent(monomial)
            shearIndeces = [3 4 5];
            itHasNoShear = ~any(monomial(shearIndeces));
        end        

        function str = num2strIncludingDots(num)
            str = strrep(num2str(num),'.','');  
        end        
        
    end    
    
end