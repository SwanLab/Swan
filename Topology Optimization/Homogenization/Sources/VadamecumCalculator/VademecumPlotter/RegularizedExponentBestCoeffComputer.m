classdef RegularizedExponentBestCoeffComputer < handle
    
   properties (Access = private)
        aValues
        bValues
        rValues
        qValues
        aValue
        bValue
        rValue
        qValue        
        nA
        nB
        nR
        nQ
        nT
        error
    end
    
    properties (Access = private)
       errorComputer 
    end
    
    methods (Access = public)
        
        function obj = RegularizedExponentBestCoeffComputer(cParams)
            obj.init(cParams);    
            obj.initValues();            
            obj.computeAxisValues();
        end
        
        function [a,b,r,q] = compute(obj)
            obj.evaluateAllSamples();
            [a,b,r,q] = obj.obtainBestParameterSample();
        end
        
        function [a,b,r,q] = obtainBestParameters(obj)
           a = 23.620689655172413;
           b = 1.444444444444444;
           r = 0.677777777777778;
           q = 11.322033898305085;              
        end          
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.errorComputer = cParams.errorComputer;
            obj.nA = 30;
            obj.nB = 10;
            obj.nR = 10;
            obj.nQ = 60; 
            obj.nT = obj.nA*obj.nB*obj.nR*obj.nQ;
        end        
        
        function initValues(obj)
            obj.aValue = zeros(obj.nT,1);
            obj.bValue = zeros(obj.nT,1);
            obj.rValue = zeros(obj.nT,1);
            obj.qValue = zeros(obj.nT,1);
            obj.error = zeros(obj.nT,1);            
        end
        
        function computeAxisValues(obj)
            obj.aValues = linspace(22,25,obj.nA);
            obj.bValues = linspace(1.2,1.6,obj.nB);
            obj.rValues = linspace(0.65,0.68,obj.nR);
            obj.qValues = linspace(11,12,obj.nQ);             
        end
        
        function evaluateAllSamples(obj)
            iter = 1;
            for iR = 1:length(obj.rValues)
                for iQ = 1:length(obj.qValues)
                    for iB = 1:length(obj.bValues)
                        for iA = 1:length(obj.aValues)     
                            obj.obtainSampleValue(iter,iA,iB,iR,iQ);
                            obj.computeError(iter);
                            iter = iter + 1;                            
                        end
                    end
                end
            end               
        end
        
        function obtainSampleValue(obj,iT,iA,iB,iR,iQ)
            obj.aValue(iT) = obj.aValues(iA);
            obj.bValue(iT) = obj.bValues(iB);
            obj.rValue(iT) = obj.rValues(iR);
            obj.qValue(iT) = obj.qValues(iQ);
        end
        
        function  computeError(obj,iT)
           a = obj.aValue(iT);
           b = obj.bValue(iT);
           r = obj.rValue(iT);
           q = obj.qValue(iT);
           e = obj.errorComputer(a,b,r,q);                            
           obj.error(iT) = e;
        end
        
        function [a,b,r,q] = obtainBestParameterSample(obj)
            [~,ind] = min(obj.error);
            a = obj.aValue(ind);
            b = obj.bValue(ind);
            r = obj.rValue(ind);
            q = obj.qValue(ind);             
        end        
               
    end
    
end