classdef testFourthOrderAmplificatorTensor < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testShapeStressWithAmplificator';
    end
    
    properties (Access = private)  
        pNorm
        monomials
        PinfRectangle
        P2Rectangle
        PinfSmoothRectangle
        P2SmoothRectangle
    end
    
    
    methods (Access = public)
        
        function obj = testFourthOrderAmplificatorTensor()
            obj.init();
            obj.computeSmoothAmplificators();
            obj.computeNonSmoothAmplificators();
            obj.plotNorms();                       
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.pNorm = 8;
        end
        
        function computeSmoothAmplificators(obj)
            d.microCase = 'SmoothRectangle';
            d.pNorm = obj.pNorm;
            aC = AmplificatorComponentsMeasurer(d);
            aC.compute();
            obj.PinfSmoothRectangle = aC.PmaxNorm;
            obj.P2SmoothRectangle   = aC.Pl2Norm;
            obj.monomials = aC.monomials;
        end
        
        function computeNonSmoothAmplificators(obj)
            d.microCase = 'Rectangle';
            d.pNorm = obj.pNorm;            
            aC = AmplificatorComponentsMeasurer(d);
            aC.compute();
            obj.PinfRectangle = aC.PmaxNorm;
            obj.P2Rectangle = aC.Pl2Norm;
            obj.monomials = aC.monomials;            
        end        
        
        function plotNorms(obj)
            obj.plotNorm([obj.PinfSmoothRectangle,obj.PinfRectangle],'PtensorCompLogInfNorm');
            obj.plotNorm([obj.P2SmoothRectangle,obj.P2Rectangle],'PtensorCompLogL2Norm');
        end
        
        function plotNorm(obj,norm,plotName)
            fig = figure();
            h = bar(log(norm));
            monomialsLeg = obj.computeAlphaLegend();
            set(gca, 'XTickLabel',monomialsLeg, 'XTick',1:numel(monomialsLeg));
            set(gca,'XTickLabelRotation',45);
            legend('SmoothRectangle','Rectangle')
            p = barPrinter(fig,h);
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            p.print([path,plotName]);
        end
        
        function a = computeAlphaLegend(obj)
            for ia = 1:size(obj.monomials,1)
               a{ia} = mat2str(obj.monomials(ia,:)); 
            end
        end
        
    end
    
    
end