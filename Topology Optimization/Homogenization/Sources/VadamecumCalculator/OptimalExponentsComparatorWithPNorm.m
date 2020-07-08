classdef OptimalExponentsComparatorWithPNorm < handle
    
    properties (Access = private)
        optimalComputer        
    end
    
    properties (Access = private)
        optimalExponentParams
        pNorm
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentsComparatorWithPNorm(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeOptimalExponents();
            obj.plotStressNormRelationWithQ();
            obj.plotConvergence();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.optimalExponentParams.txi  = pi/4;
            obj.optimalExponentParams.rho  = 0.15;
            obj.optimalExponentParams.phi  = pi/4;
            obj.optimalExponentParams.hMesh = 0.01;
            obj.pNorm = {2,4,16,32,'max'};            
        end
        
        function computeOptimalExponents(obj)            
            s = obj.optimalExponentParams;            
            for iP = 1:length(obj.pNorm)
                p = obj.pNorm{iP};
                s.pNorm = p;
                s.fileName = ['ExamplePaper',num2str(p)];    
                oC = OneOptimalExponentComputerAndFunctionVariation(s);
                oC.compute();
                obj.optimalComputer{iP} = oC;
            end
            
        end
        
        function plotStressNormRelationWithQ(obj)
            f = figure();
            hold on
            for iP = 1:length(obj.pNorm)
                h{iP} = plot(obj.optimalComputer{iP}.qValues(2:end),obj.optimalComputer{iP}.fValues(2:end),'-+');
                p = obj.pNorm{iP};                
                if isequal(p,'max')
                    leg{iP} = '$p = \infty$';
                else
                    leg{iP} = ['p = ',num2str(p)];
                end                
            end
            outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
            fTitle = 'Stress norm';
            xlabel('$q$','interpreter','latex')
            % set(gca, 'XTickLabel',[])
            ylabel(fTitle,'interpreter','latex')
            
            for iP = 1:length(obj.pNorm)
                x = obj.optimalComputer{iP}.qValues(2:end);
                y = obj.optimalComputer{iP}.fValues(2:end);                
                [~,imin] = min(y);
                plot(x(imin),y(imin),'-s','Color',h{iP}.Color,'MarkerSize',10,'MarkerEdgeColor',h{iP}.Color,'MarkerFaceColor',h{iP}.Color)
            end            
            
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');             
            %outputName = [outPutPath,'StressNormVsQForPmax'];
            %outputName = [outPutPath,'StressNormVsQForDifferentP'];
            printer = plotPrinter(f,h);
           % printer.print(outputName);
        end
        
        function plotConvergence(obj)
            f = figure();
            hold on
            for iP = 1:length(obj.pNorm)           
                h{iP} = plot(obj.optimalComputer{iP}.fOptIter,'-+');
                p = obj.pNorm{iP};
                if isequal(p,'max')
                    leg{iP} = '$p = \infty$';
                else
                    leg{iP} = ['p = ',num2str(p)];
                end
            end
            outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
            fTitle = 'Stress norm';
            xlabel('iterations','interpreter','latex')
            ylabel(fTitle,'interpreter','latex')
            
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');                       
            outputName = [outPutPath,'StressNormMinimizationForPmax'];
            %outputName = [outPutPath,'StressNormMinimizationFineMax'];
            printer = plotPrinter(f,h);
           % printer.print(outputName);
        end
        
        
        
    end
    
    
    
    
    
end