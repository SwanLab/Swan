classdef OptimalExponentComputer < handle
    
    properties (Access = public)
       qOpt 
       qMin
       qMax
    end
    
    properties (Access = private)
        maxStress
        mx
        my
        rho
        txi
        rhoV
        txiV
        mxMax
        myMax
        phi
        phiV
        samplePoints
        fileName
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            nphi = length(obj.phiV);
            npoints = length(obj.rhoV);
            obj.qOpt = zeros(nphi,npoints);
            obj.qMin = zeros(nphi,npoints);
            obj.qMax = zeros(nphi,npoints);
            for iphi = 1:nphi
                obj.phi = obj.phiV(iphi);
                disp([num2str(iphi/nphi*100),'%'])                
                for ipoint = 1:npoints
                    obj.rho = obj.rhoV(ipoint);
                    obj.txi = obj.txiV(ipoint);
                    [q,qmin,qmax] = obj.computeOptimalExponent(ipoint,iphi);
                    obj.qOpt(iphi,ipoint) = q;
                    obj.qMin(iphi,ipoint) = qmin;
                    obj.qMax(iphi,ipoint) = qmax;
                end
            end
            x.rho = obj.rhoV;
            x.txi = obj.txiV;
            x.phi = obj.phiV;
            x.q = obj.qOpt;
            save(obj.fileName,'x');
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.samplePoints = cParams.samplePoints;
            obj.computeRhoTxiAndPsiSamples();
        end
        
        function computeRhoTxiAndPsiSamples(obj)
            obj.rhoV = obj.samplePoints.rhoV;
            obj.txiV = obj.samplePoints.txiV;
            obj.phiV = obj.samplePoints.phiV;
        end
        
        function [q,qMin,qMax] = computeOptimalExponent(obj,ipoint,iphi)            
            s.fileName = [obj.fileName,'Txi',num2str(ipoint),'Phi',num2str(iphi)];
            s.rho   = obj.rho;
            s.txi   = obj.txi;
            s.phi   = obj.phi;
            s.pNorm = 'max';
            s.hMesh = 0.1;
            c = OneOptimalExponentComputerAndFunctionVariation(s);
            c.computeOptimalExponent();
            %c.printOptimalMicroStructure();
            
            qMin = c.qOptIter();
            fMin = c.fOptIter();
            
            [~,ind] = min(fMin);
            q = qMin(ind);
            qMin = c.qMax;
            qMax = c.qMin;
        end

    end    
    
end



