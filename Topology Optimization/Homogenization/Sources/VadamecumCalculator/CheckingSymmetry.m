classdef CheckingSymmetry < handle
    
    properties (Access = private)
       phi
       rho
       txi
       qValues
       fValues
       stressNorm
       symmetricError
    end
    
    properties (Access = private)
       phiV 
       txiV
       rhoV      
       fileName
       outputPath
       nQ
    end
    
    methods (Access = public)
        
        
        function obj = CheckingSymmetry()
            obj.init();
            obj.compute();
        end
 
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.phiV = linspace(0,pi,30);
            obj.rhoV  = [0.9,0.9,0.5,0.5];
            obj.txiV  = pi/2 - [0.1083,0.557,0.88974,1.0984];
            obj.fileName = 'CheckingSymmetry';
            obj.outputPath = '/home/alex/git-repos/MicroStructurePaper/';            
            obj.nQ = 10;
        end
        
        function compute(obj)
            for iTest = 1:length(obj.rhoV)
                obj.rho   = obj.rhoV(iTest);
                obj.txi   = obj.txiV(iTest);
                for iphi = 1:length(obj.phiV)
                    obj.phi = obj.phiV(iTest);
                    obj.computeSymetricCases();
                    obj.computeSymmetryError(iTest);
                    obj.plotSymmetryError(iTest,iphi);
                end
            end
        end
        
        function computeSymetricCases(obj)    
            [q1,f1] = obj.computeOptimalValue(obj.txi,obj.rho,obj.phi);
            [q2,f2] = obj.computeOptimalValue(pi/2 - obj.txi,obj.rho,pi/2 -obj.phi);
            obj.qValues = [q1 q2];
            obj.fValues = [f1 f2];
        end
        
        function [qV,fV] = computeOptimalValue(obj,txi,rho,phi)
            phi = obj.shiftPhiInInterval(phi);
            obj.createStressNormCalculator(txi,rho,phi);
            qV(:,1) = linspace(obj.stressNorm.qMin,obj.stressNorm.qMax,obj.nQ);                        
            fV = zeros(obj.nQ,1);
            for iq = 1:obj.nQ
                fV(iq,1) = obj.stressNorm.objective(qV(iq));
            end
        end
        
        function createStressNormCalculator(obj,txi,rho,phi)
            s.fileName = obj.fileName;
            s.rho   = rho;
            s.xi   = txi;
            s.phi   = phi;
            s.pNorm = 16;
            s.hMesh = 0.1;
            s.print = false;
            s.hasToCaptureImage = false;
            c = StressNormVsQproblemCreator(s);
            c.compute();            
            obj.stressNorm = c;
        end
        
        function computeSymmetryError(obj,iTest)
            q = obj.qValues;
            f = obj.fValues;
            errorQ = norm(q(:,1) - q(:,2))/norm(q(:,1));
            errorF = norm(f(:,1) - f(:,2))/norm(f(:,1));
            obj.symmetricError(iTest) = errorQ + errorF;
        end
        
        function plotSymmetryError(obj,iTest,iPhi)            
            f = figure();
            hold on;
            h{1} = plot(obj.qValues(:,1),obj.fValues(:,1),'-+','LineWidth',4);
            h{2} = plot(obj.qValues(:,2),obj.fValues(:,2),'-+','LineWidth',4);
            xlabel('$q$','Interpreter','latex');
            ylabel('$||\sigma||^p$','Interpreter','latex');
            obj.addLegend();
            obj.addTitle();
            obj.print(f,iTest,iPhi);
        end
        
        function addLegend(obj)
            l1 = '$||\sigma||^p(q,\xi,\phi)$';
            l2 = '$||\sigma||^p(q,\pi - \xi,\pi -\phi)$'; 
            legend(l1,l2,'Interpreter','latex','Location','best');                        
        end
        
        function addTitle(obj)
            rhoT = ['\rho = ',num2str(obj.rho)];
            txiT = ['\xi = ',num2str(obj.txi)];
            phiT = ['\phi = ',num2str(obj.phi)];            
            tit =  ['$',rhoT,'\quad',txiT,'\quad',phiT,'$'];
            title(tit,'Interpreter','latex');            
        end
        
        function print(obj,f,iTest,iPhi)
            fp = contourPrinter(f);
            filePath = [obj.outputPath,'CheckingSymmetry',num2str(iTest),num2str(iPhi)];
            fp.print(filePath);            
        end
        
        function phi = shiftPhiInInterval(obj,phi)            
            phi = wrapTo2Pi(phi); 
            if phi > pi
               phi = wrapTo2Pi(phi + pi); 
            end            
        end        
        
    end
        
    
end