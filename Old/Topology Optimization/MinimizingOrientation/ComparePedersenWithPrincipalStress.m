classdef ComparePedersenWithPrincipalStress < handle
    
    properties (Access = private)
        Ctensor
        stress
        principalDirection
        optimalAngle
        minimalCompOrientation
        gamma
        coef
        cases
        icase
        orientations
        principalStress
        alphas
        eta
        fig
        hFig
        pointWidth
        constitutiveType
        stressCases
        tensorCases
        iStress
        iTensor
    end
    
    methods (Access = public)
        
        function obj = ComparePedersenWithPrincipalStress()
            obj.init();
            for iStressCase = 1:numel(obj.stressCases)
                for iTensorCase = 1:numel(obj.tensorCases)
                    obj.iStress = iStressCase;
                    obj.iTensor = iTensorCase;
                    obj.createStress();
                    obj.createConstitutiveTensor();
                    obj.computeConstitutiveCoeficient();
                    obj.computePrincipalStress();
                    obj.computeEta();
                    obj.computeGamma();
                    obj.computeOrientations();
                    obj.computeMinimalComplianceOrientation();
                    obj.compareResults();
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.cases = {'MinimumEigenValue',...
                'MaximumEigenValue',...
                'MinimumAbsEigenValue',...
                'MaximumAbsEigenValue'};
            obj.stressCases = {'First','Second','Third'};
            obj.tensorCases = {'Vertical','Horizontal','Square','Isotropic'};
        end
        
        function createConstitutiveTensor(obj)
            switch obj.tensorCases{obj.iTensor}
                case 'RandomOrthotropic'
                    obj.createOrthotropicConstitutiveTensor();
                case 'Vertical'
                    obj.createVerticalOrthotropicConstitutiveTensor();
                case 'Horizontal'
                    obj.createHorizontalOrthotropicConstitutiveTensor();
                case 'Square'
                    obj.createSquareOrthotropicConstitutiveTensor();
                case 'Isotropic'
                    obj.createIsotropicConstitutiveTensor();
            end
            
        end
        
        function createStress(obj)
            switch obj.stressCases{obj.iStress}
                case 'First'
                    obj.stress = [1 -2 0.5]';
                case 'Second'
                    obj.stress = [0.7,0.7,0.3]';
                case 'Third'
                    obj.stress = [1 -2 0.0001]';
            end
        end
        
        function createOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = rand(1,1);
            C(2,2) = rand(1,1);
            C(1,2) = obj.randBetweenBounds(-0.1,0.1);
            C(2,1) = C(1,2);
            C(3,3) = obj.randBetweenBounds(0,1);
            obj.Ctensor = C;
        end
        
        function createVerticalOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = 0.374827;
            C(2,2) = 0.864285;
            C(1,2) = 0.1257865;
            C(2,1) = C(1,2);
            C(3,3) = 0.13265;
            obj.Ctensor = C;
        end
        
        function createHorizontalOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = 0.864285;
            C(2,2) = 0.374827;
            C(1,2) = 0.1257865;
            C(2,1) = C(1,2);
            C(3,3) = 0.13265;
            obj.Ctensor = C;
        end
        
        function createSquareOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = 0.664285;
            C(2,2) = 0.664285;
            C(1,2) = 0.1257865;
            C(2,1) = C(1,2);
            C(3,3) = 0.13265;
            obj.Ctensor = C;
        end
        
        function createIsotropicConstitutiveTensor(obj)
            E = obj.randBetweenBounds(0.01,1);
            nu = obj.randBetweenBounds(0,0.49);
            C = zeros(3,3);
            C(1,1) = 1;
            C(2,2) = 1;
            C(1,2) = nu;
            C(2,1) = C(1,2);
            C(3,3) = (1-nu)/2;
            C = E*C;
            obj.Ctensor = C;
        end
        
        function computeConstitutiveCoeficient(obj)
            C = obj.Ctensor;
            C11 = C(1,1);
            C22 = C(2,2);
            C12 = C(1,2);
            C66 = C(3,3);
            obj.coef = (C11 - C22)/((C11+C22)-2*(C12 + 2*C66));
        end
        
        function computePrincipalStress(obj)
            cParams.type = '2D';
            p = PrincipalDirectionComputer.create(cParams);
            s(1,:) = obj.stress;
            p.compute(s);
            obj.principalStress = p;
        end
        
        function computeOrientations(obj)
            for ic = 1:numel(obj.cases)
                obj.icase = ic;
                obj.computeOrientationCase();
            end
        end
        
        
        function computeOrientationCase(obj,cParams)
            cParams.type = obj.cases{obj.icase};
            oUpdater = OrientationUpdater.create(cParams);
            cParams.pD = obj.principalStress.direction;
            cParams.pS = obj.principalStress.principalStress;
            oUpdater.compute(cParams);
            obj.orientations{obj.icase} = oUpdater.alpha;
            obj.alphas{obj.icase} = atan2(oUpdater.alpha(2),oUpdater.alpha(1));
            %obj.alphas{obj.icase} = acos(oUpdater.alpha(1))*sign(oUpdater.alpha(2));
        end
        
        function computeMinimalComplianceOrientation(obj)
            p = OptimalOrientationComputer();
            p.compute(obj.stress,obj.Ctensor);
            obj.minimalCompOrientation = p;
            obj.optimalAngle = p.optimalAngle;
        end
        
        function plotAllDomain(obj)
            comp = obj.minimalCompOrientation.compliance;
            xp = linspace(-pi,pi,100);
            obj.hFig{1} = plot(xp*180/pi,comp(xp),'--');
            set(obj.hFig{end},'LineWidth',1.5);
        end
        
        function plotExtremeValues(obj)
            comp = obj.minimalCompOrientation.compliance;
            for ic = 1:numel(obj.cases)
                x = obj.alphas{ic};
                obj.hFig{end+1} = plot(x*180/pi,comp(x),'+');
                set(obj.hFig{end},'LineWidth',obj.pointWidth);
            end
        end
        
        function plotOptimalValue(obj)
            comp = obj.minimalCompOrientation.compliance;
            x = obj.optimalAngle;
            obj.hFig{end+1} = plot(x*180/pi,comp(x),'+');
            set(obj.hFig{end},'LineWidth',obj.pointWidth);
        end
        
        function compareResults(obj)
            obj.pointWidth = 4;
            obj.fig = figure();
            hold on
            obj.plotAllDomain();
            obj.plotExtremeValues();
            %obj.plotOptimalValue();
            obj.addNames();
            obj.addLegend();
            obj.addTitle();
            obj.setPositionAndSize();
            obj.printFigure();
        end
        
        function setPositionAndSize(obj)
            set(obj.fig,'units','points','position',[1920,500,600,600])
        end
        
        function addNames(obj)
            xlabel('\alpha')
            ca = get(obj.fig,'CurrentAxes');
            xl = get(ca,'xlabel');
            set(xl,'FontSize',30)
        end
        
        function addTitle(obj)
            stressStr = [];
            for i = 1:length(obj.stress)
                str = num2str(obj.stress(i));
                nStr = [str(~isspace(str)),', '];
                stressStr = [stressStr,nStr];
            end
            t = title([obj.tensorCases{obj.iTensor}, ' with \sigma = [',stressStr(1:end-2),']']);
            t.FontSize = 20;
        end
        
        function addLegend(obj)
            legendName{1} = 'All' ;
            for ic = 1:numel(obj.cases)
                legendName{end+1} = obj.cases{ic};
            end
            %legendName{end+1} = 'Optimal';
            %legend(legendName,'Location','Best');
            legend(legendName)
            legend(legendName,'Location','North');
        end
        
        function printFigure(obj)
            path = '/home/alex/Dropbox/Amplificators/GregoireMeeting7/OptimalOrientation';
            figName = [path,obj.tensorCases{obj.iTensor},obj.stressCases{obj.iStress}];
            print(obj.fig,figName,'-dpng');
        end
        
        function computeEta(obj)
            stressPrin = obj.principalStress.principalStress;
            s = stressPrin;
            if abs(s(1)) > abs(s(2))
                t = s(2)/s(1);
            else
                t = s(1)/s(2);
            end
            obj.eta = t;
        end
        
        function computeGamma(obj)
            obj.gamma = obj.coef*obj.eta;
        end
        
    end
    
    methods (Access = private, Static)
        
        function c = randBetweenBounds(a,b)
            c = (b-a).*rand(1,1) + a;
        end
        
    end
    
    
    
end