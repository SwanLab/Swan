classdef OptimalExponentsComparatorWithHmesh < handle
    
    properties (Access = private)
        optimalComputer        
    end
    
    properties (Access = private)
        optimalExponentParams
        hMesh
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentsComparatorWithHmesh(cParams)
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
            obj.hMesh = [0.01,0.005,0.0025];
            obj.optimalExponentParams.pNorm = 'max';            
        end
        
        function computeOptimalExponents(obj)            
            s = obj.optimalExponentParams;            
            for iH = 1:length(obj.hMesh)
                h = obj.hMesh(iH);
                s.hMesh = h;
                s.fileName = ['ExamplePaper',strrep(num2str(h),'.','_')];    
                oC = OneOptimalExponentComputerAndFunctionVariation(s);
                oC.compute();
                obj.optimalComputer{iH} = oC;
            end
            
        end
        
        function plotStressNormRelationWithQ(obj)
            f = figure();
            hold on
            for ih = 1:length(obj.hMesh)
                h{ih} = plot(obj.optimalComputer{ih}.qValues(2:end),obj.optimalComputer{ih}.fValues(2:end),'-+');
                hM = obj.hMesh(ih);                
                leg{ih} = ['h = ',num2str(hM)];
            end
            outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
            fTitle = 'Stress norm';
            xlabel('$q$','interpreter','latex')
            % set(gca, 'XTickLabel',[])
            ylabel(fTitle,'interpreter','latex')
           
            for ih = 1:length(obj.hMesh)
                x = obj.optimalComputer{ih}.qValues(2:end);
                y = obj.optimalComputer{ih}.fValues(2:end);                
                [~,imin] = min(y);
                plot(x(imin),y(imin),'-s','Color',h{ih}.Color,'MarkerSize',10,'MarkerEdgeColor',h{ih}.Color,'MarkerFaceColor',h{ih}.Color)
            end
            
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');             
            %outputName = [outPutPath,'StressNormVsQ'];
                      
            outputName = [outPutPath,'StressNormVsQWithH'];
            printer = plotPrinter(f,h);
           % printer.print(outputName);
        end
        
        function plotConvergence(obj)
            f = figure();
            hold on
            for iH = 1:length(obj.hMesh)           
                h{iH} = plot(obj.optimalComputer{iH}.fOptIter,'-+');
                hM = obj.hMesh(iH);
                leg{iH} = ['h = ',num2str(hM)];
            end
            outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
            fTitle = 'Stress norm';
            xlabel('iterations','interpreter','latex')
            ylabel(fTitle,'interpreter','latex')
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');                       
            %outputName = [outPutPath,'StressNormMinimization'];
            outputName = [outPutPath,'StressNormMinimizationMaxWithH'];
            printer = plotPrinter(f,h);
            %printer.print(outputName);
        end
        
    end
    
    
end