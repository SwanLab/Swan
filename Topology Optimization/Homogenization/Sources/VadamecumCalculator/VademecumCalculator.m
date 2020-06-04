classdef VademecumCalculator < handle
    
    properties (Access = private)
        samplePoints
        optimalExponent
        iter
        iMx
        iMy   
        cellVariables
        qOptimal
        phi
    end
    
    properties (Access = private)
        nMx
        nMy
        nPhi
        phiMin
        phiMax
    end
    
    methods (Access = public)
        
        function obj = VademecumCalculator()
            obj.init([])
            obj.createSamplePoints();
            obj.computeVademecumData();
        end
        
        function computeVademecumData(obj)
            for imx = 1:obj.nMx
                for imy = 1:obj.nMy
                    obj.storeIndeces(imx,imy);
                    obj.computeGlobalIteration();
                    obj.displayPercentatge();
                    for iphi = 1:obj.nPhi                               
                        obj.createOptimalExponentComputer(iphi);
                        obj.computeOptimalExponent();
                        obj.obtainCellVariables(iphi);
                    end
                end
            end
        end
        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nMx    = 10;
            obj.nMy    = 10;
            obj.nPhi   = 10;
            obj.phiMin = 0;
            obj.phiMax = pi/4;
        end
        
        function createSamplePoints(obj)
            s.type = 'FromMxMy';
            s.nMx  = obj.nMx;
            s.nMy  = obj.nMy;
            s.nPhi = obj.nPhi;
            s.phiMin = obj.phiMin;
            s.phiMax = obj.phiMax;
            sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            sample.compute();
            obj.samplePoints = sample;
        end
        
        function storeIndeces(obj,imx,imy)
            obj.iMx = imx;
            obj.iMy = imy;
        end
        
        function computeGlobalIteration(obj)
            imx = obj.iMx;
            imy = obj.iMy;
            obj.iter = (imy + obj.nMx*(imx -1));
        end

        function displayPercentatge(obj)
            disp([num2str(obj.iter/(obj.nMx*obj.nMy)*100),'% done']);
        end       
        
        function q = computeOptimalExponent(obj)
             obj.optimalExponent.computeOptimalExponent();  
             qMin = obj.optimalExponent.qOptIter();
             fMin = obj.optimalExponent.fOptIter();            
             [~,ind] = min(fMin);
             q = qMin(ind);
             %q = 2;
             obj.qOptimal(obj.iMx,obj.iMy) = q;
        end
        
        function obtainCellVariables(obj,iphi)
            imx = obj.iMx;
            imy = obj.iMy;
            q = obj.qOptimal(imx,imy);
            cVariables = obj.optimalExponent.obtainCellVariables(q);  
            cVariables.rho = obj.samplePoints.rhoV(obj.iter);
            cVariables.xi = obj.samplePoints.txiV(obj.iter);
            cVariables.phi = obj.samplePoints.phiV(iphi);
            cVariables.q = q;
            obj.cellVariables{imx,imy} = cVariables;
        end
        
        function createOptimalExponentComputer(obj,iphi)
            s.fileName = ['OptimaSuperEllipseIter',num2str(obj.iter),'Phi',num2str(iphi)];
            s.rho   = obj.samplePoints.rhoV(obj.iter);
            s.txi   = obj.samplePoints.txiV(obj.iter);
            s.phi   = obj.samplePoints.phiV(iphi);
            s.pNorm = 'max';
            s.hMesh = [];            
            optimalExp = OneOptimalExponentComputerAndFunctionVariation(s);           
            obj.optimalExponent = optimalExp;            
        end
        
    end
    
    
    
end