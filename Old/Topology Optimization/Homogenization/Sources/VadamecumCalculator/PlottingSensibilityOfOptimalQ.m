classdef PlottingSensibilityOfOptimalQ < handle
    
    properties (Access = private)
        stressProblem
        rho
        xi
        phi
        fValues
        qValues
        qOpt
        fOpt
        hPlot
        figureN
    end
    
    properties (Access = private)
        rhoV
        xiV
        phiV
        fileName
        outPutPath
        hMesh
        pNorm
        nValues
    end
    
    
    methods (Access = public)
        
        function obj = PlottingSensibilityOfOptimalQ()
            obj.init();
            for iTest = 3:length(obj.rhoV)
                obj.rho = obj.rhoV(iTest);
                obj.xi  = obj.xiV(iTest);
                obj.computePhiV();
                for iphi = 1:length(obj.phiV)
                    obj.phi = obj.phiV(iphi);
                    obj.createStressProblem();
                    obj.computeStressNormsVsQ(iphi);
                    obj.computeOptimalQ(iphi);
                end
                obj.plotImage();
                obj.printImage(iTest);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.rhoV = [0.9 0.9 0.5 0.5];
            obj.xiV  = [83.7951 58.0865 39.0219 27.0665]*pi/180;
            obj.hMesh = 0.2;
            obj.pNorm = 'max';
            obj.fileName = 'MaxStressVsQ';
            obj.outPutPath = '/home/alex/Dropbox/GregMeeting/';
            obj.nValues = 50;
        end
        
        function computePhiV(obj)
            %obj.phiV = [20*pi/32,21*pi/32,22*pi/32,23*pi/32,24*pi/32,...
            %            obj.xi,pi - obj.xi,0,pi/4,pi/2,3*pi/4];       
            obj.phiV = [obj.xi,pi - obj.xi,0,pi/4,pi/2,3*pi/4];                      
        end
        
        function createStressProblem(obj)
            s.rho = obj.rho;
            s.txi = obj.xi;
            s.fileName = obj.fileName;
            s.phi = obj.phi;
            s.hMesh = obj.hMesh;
            s.pNorm = obj.pNorm;
            sProblem = OneOptimalExponentComputerAndFunctionVariation(s);
            obj.stressProblem = sProblem;
        end
        
        function computeStressNormsVsQ(obj,iter)
            obj.stressProblem.computeStressNormRelationWithQ(obj.nValues);
            obj.fValues{iter} = obj.stressProblem.fValues;
            obj.qValues{iter} = obj.stressProblem.qValues;
        end
        
        function computeOptimalQ(obj,iter)
            obj.stressProblem.computeOptimalExponent();
            obj.qOpt{iter} = obj.stressProblem.qOptIter;
            obj.fOpt{iter} = obj.stressProblem.fOptIter;
        end
        
        function plotImage(obj)
            obj.figureN = figure();
            hold on
            obj.plotStressNormsVsQ();
            obj.plotOptimalQ();
            obj.addAxisLimits();
            obj.addLegend();
            obj.addTitle();
        end
        
        function h = plotStressNormsVsQ(obj)
            nPhi = length(obj.phiV);
            h = cell(nPhi,1);
            for iphi = 1:nPhi
                q = obj.qValues{iphi};
                f = obj.fValues{iphi};
                h{iphi} = plot(q(1:end),f(1:end),'-+','LineWidth',1);
            end
            obj.hPlot = h;
        end
        
        function plotOptimalQ(obj)
            h = obj.hPlot;
            for iphi = 1:length(obj.phiV)
                q = obj.qOpt{iphi};
                f = obj.fOpt{iphi};
                [~,ind] = min(f);
                plot(q(ind),f(ind),'-s','LineWidth',2,'Color',h{iphi}.Color)
            end
        end
        
        function addAxisLimits(obj)
            ax = gca;
            ax.XLim = [2 32];
        end
        
        function addLegend(obj)
            l1 = '$\phi = \xi$';
            l2 = '$\phi = \pi -\xi$';
            l3 = '$\phi = 0$';
            l4 = '$\phi = \pi/4$';
            l5 = '$\phi = \pi/2$';
            l6 = '$\phi = 3\pi/4$';
            legend(l1,l2,l3,l4,l5,l6,'Interpreter','latex','Location','northeast');
        end
        
        function addTitle(obj)
            xiT  = ['\xi = ',num2str(obj.xi)];            
            rhoT = ['\rho = ',num2str(obj.rho)];            
            title(['$',rhoT,',\quad ',xiT,'$'],'Interpreter','latex');
        end
        
        function printImage(obj,iTest)
            f = obj.figureN;
            h = obj.hPlot;
            outputName = [obj.outPutPath,obj.fileName,num2str(iTest)];
            printer = plotPrinter(f,h);
            printer.print(outputName);
        end
        
    end
    
    
end