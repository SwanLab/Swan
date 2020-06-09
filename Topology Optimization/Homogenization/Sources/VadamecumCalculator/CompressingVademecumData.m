classdef CompressingVademecumData < handle
    
    properties (Access = private)
        qV
        rhoV
        xiV
        phiV
        cellDataSameSign
        cellDataDifferentSign
        phiIndeces
    end
    
    properties (Access = private)
        compressedFileName
        nMx
        nMy
        nPhiQuarter
        nPhiT
    end
    
    methods (Access = public)
        
        function obj = CompressingVademecumData()
            obj.init();
            obj.initVariables();
            obj.loadVademecumDataBase();    
            obj.compute();
            obj.saveData();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.compressedFileName = 'OptimalSuperEllipses';
            obj.nMx = 20;
            obj.nMy = 20;
            obj.nPhiQuarter = 10;
            obj.nPhiT = obj.computeNphiT();
        end
        
        function nPT = computeNphiT(obj)
            n0  = obj.nPhiQuarter;
            n1  = 2*(n0 -1)+1;
            nPT = 2*(n1 -1)+1;
        end
        
        function initVariables(obj)
            nmx = obj.nMx;
            nmy = obj.nMy;
            np = obj.nPhiT;
            obj.qV   = zeros(nmx,nmy,np);
            obj.rhoV = zeros(nmx,nmy,np);
            obj.xiV  = zeros(nmx,nmy,np);
            obj.phiV = zeros(nmx,nmy,np);
        end
      
        function loadVademecumDataBase(obj)
            s.nMx = obj.nMx;
            s.nMy = obj.nMy;
            s.nPhi = obj.nPhiQuarter;
            v = VademecumDataLoader(s);
            obj.cellDataSameSign = v.cellDataSameSign;
            obj.cellDataDifferentSign = v.cellDataDifferentSign;
        end
        
        function compute(obj)            
            for iMx = 1:obj.nMx
                for iMy = 1:obj.nMy
                    for iPhi = 1:nPhi
                        obj.computePhiIndices(iPhi)
                        obj.computeXiValues(iMx,iMy,iPhi);
                        obj.computeRhoValues(iMx,iMy,iPhi);
                        obj.computePhiValues(iMx,iMy,iPhi);
                        obj.computeQValues(iMx,iMy,iPhi);
                        obj.plotFigures(iMx,iMy,iPhi);
                    end
                end
            end           
        end
        
        function computePhiIndices(obj,iPhi)            
            iphi1 = iPhi;
            iphi2 = 2*nPhi - iPhi;
            iphi3 = 2*(nPhi -1) + iPhi;
            iphi4 = 2*(2*nPhi -1) - iPhi;            
            obj.phiIndeces = [iphi1,iphi2,iphi3,iphi4];
        end
        
        function computeXiValues(obj,iMx,iMy,iPhi)
             [cS,cD] = obj.obtainCellData(obj,iMx,iMy,iPhi);
             xiQ = [cS.xi cS.xi cD.xi cD.xi];
             ind = obj.phiIndeces;
             obj.xiV(iMx,iMy,ind) = xiQ;
        end
        
        function computeRhoValues(obj,iMx,iMy,iPhi)
             [cS,cD] = obj.obtainCellData(obj,iMx,iMy,iPhi);
             rhoQ = [cS.rho cS.rho cD.rho cD.rho];
             ind = obj.phiIndeces;
             obj.rhoV(iMx,iMy,ind) = rhoQ;
        end
        
        function computePhiValues(obj,iMx,iMy,iPhi)
            [cS,cD] = obj.obtainCellData(obj,iMx,iMy,iPhi);         
            qQ = [cS.phi pi/2-cS.phi cD.phi 3*pi/2-cD.phi];
            ind = obj.phiIndeces;
            obj.qV(iMx,iMy,ind) = qQ;
        end                
        
        function computeQValues(obj,iMx,iMy,iPhi)
            [cS1,cD1] = obj.obtainCellData(obj,iMx,iMy,iPhi);
            [cS2,cD2] = obj.obtainCellData(obj,iMy,iMx,iPhi);
            qQ = [cS1.q cS2.q cD1.q cD2.q];
            ind = obj.phiIndeces;
            obj.qV(iMx,iMy,ind) = qQ;
        end
        
        function [cS,cD] = obtainCellData(obj,iMx,iMy,iPhi)
           cS = obj.cellDataSameSign{iMx}{iMx,iMy,iPhi};
           cD = obj.cellDataDifferentSign{iMx}{iMx,iMy,iPhi};            
        end
        
        function plotFigures(obj,iMx,iMy,iPhi)
            [cS1,cD1] = obj.obtainCellData(obj,iMx,iMy,iPhi);
            [cS2,cD2] = obj.obtainCellData(obj,iMy,iMx,iPhi);
            obj.plotFigure(cS1,1);
            obj.plotFigure(cD1,2);
            obj.plotFigure(cS2,3);
            obj.plotFigure(cD2,4);
        end
        
       function saveData(obj)
            q = obj.qV;
            rho = obj.rhoV;
            xi  = obj.xiV;
            phi = obj.phiV;
            save(obj.compressedFileName,'q','rho','xi','phi')
        end
 
    end
    
    methods (Access = private, Static)
        
        function [xiP,rhoP,phiP,qP] = loadData(cP)
            xiP  = cP.xi;
            rhoP = cP.rho;
            phiP = cP.phi;
            qP   = cP.q;          
        end
        
        function plotFigure(cellData,nFigure)
            figure(nFigure);
            clf
            cellData.mesh.plot;
            drawnow
        end        
        
    end
    
end