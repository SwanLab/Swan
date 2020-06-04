classdef PrintMicroStructuresForGaussianMean < handle
    
    properties (Access = private)
       rho
       xi
       phi
       stressProblem
    end
    
    properties (Access = private)
       rhoV
       xiV
       phiV
       hMesh
       pNorm
       fileName
       outPutPath        
    end
    
    methods (Access = public)
        
        function obj = PrintMicroStructuresForGaussianMean()
            obj.init();
            for iTest = 3:3%length(obj.rhoV)
                obj.rho = obj.rhoV(iTest);
                obj.xi = obj.xiV(iTest);
                obj.computePhiV();                
                for iphi = 1:length(obj.phiV)
                    obj.phi = obj.phiV(iphi);
                    obj.printOptimalMicroStructure(iTest,iphi);                
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.rhoV = [0.9 0.9 0.5 0.5];
            obj.xiV  = [83.7951 58.0865 39.0219 27.0665]*pi/180;
            obj.pNorm = 16;
            obj.fileName = 'OptimalMicroForStress';
            obj.outPutPath = '/home/alex/git-repos/MicroStructurePaper/';                        
        end
        
        function computePhiV(obj)
            obj.phiV = [0,pi/4,pi/2,3*pi/4,pi,obj.xi,pi - obj.xi];
        end
        
        function printOptimalMicroStructure(obj,iTest,iPhi)
            obj.createStressProblem(iTest,iPhi);
            obj.stressProblem.computeOptimalExponent();
            obj.stressProblem.printOptimalMicroStructure();
        end
        
        function createStressProblem(obj,iTest,iPhi)
            s.rho = obj.rho;
            s.txi = obj.xi;
            s.fileName = [obj.fileName,'Case',num2str(iTest),'Phi',num2str(iPhi)];
            s.phi = obj.phi;
            s.hMesh = [];
            s.pNorm = obj.pNorm;
            sProblem = OneOptimalExponentComputerAndFunctionVariation(s);
            obj.stressProblem = sProblem;
        end        
 
    end
    
end