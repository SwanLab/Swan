classdef OptimalExponentForMachingVigdergauz < handle
    
    properties (Access = private)
        sample
        rho
        txi
        rhoV
        txiV
        errorV
        qV
        mxV
        myV
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentForMachingVigdergauz()
            obj.init()
            obj.compute()
            obj.save();
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createSample();
        end
        
        function createSample(obj)
            s.type = 'FromMxMy';
            obj.sample = SamplePointsCreatorForOptimalExponentComputer.create(s);
            obj.sample.compute();
        end
        
        function compute(obj)
            for ipoint = 1: length(obj.sample.rhoV)
                obj.rho = obj.sample.rhoV(ipoint);
                obj.txi = obj.sample.txiV(ipoint);
                [qOpt,error] = obj.computeOptimalExponent();
                obj.qV(ipoint) = qOpt;
                obj.errorV(ipoint) = error;
                obj.mxV(ipoint) = obj.computeMx(obj.txi,obj.rho,qOpt);
                obj.myV(ipoint) = obj.computeMy(obj.txi,obj.rho,qOpt);                
            end            
        end
        
        function [qOpt,error] = computeOptimalExponent(obj)
            s.rho = obj.rho;
            s.txi = obj.txi;
            s.savingFrames = true;
            optimalExponent = OptimalExponentClosestToVigergauz(s);
            optimalExponent.compute();
            qOpt = optimalExponent.qOpt;
            error = optimalExponent.error;
        end
        
        function save(obj)
           d.mx = obj.mxV;
           d.my = obj.myV;
           d.rho = obj.rhoV;
           d.txi = obj.txiV;
           d.q  = obj.qV;
           d.error = obj.errorV;
           fN = 'OptimalSuperEllipseExponentMatchingVigdergauz';
           pD = 'Topology Optimization/Vademecums';
           file2SaveName = [pD,'/',fN,'.mat'];
           save(file2SaveName,'d');                    
        end

    end
    
    methods (Access = private, Static)
        
        function mx = computeMx(txi,rho,q)
            mx = SuperEllipseParamsRelator.mx(txi,rho,q);
        end
        
        function my = computeMy(txi,rho,q)
            my = SuperEllipseParamsRelator.my(txi,rho,q);            
        end
        
    end
    
    
end