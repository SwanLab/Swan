classdef AxAyComputerFromVolumeAndStrain < handle
    
    properties (Access = private)
        theta
        strain
        cx
        cy
        mu0
        mu1
        k0
        k1
    end
    
    methods (Access = public)
        
        function obj = AxAyComputerFromVolumeAndStrain(cParams)
            obj.init(cParams)
        end
        
        function [ax,ay] = compute(obj)
            b = obj.computeB();
            c = obj.computeC();
            q = b/c;
            q2 = ((obj.mu1 + obj.k0)/(obj.k1 - obj.k0) + obj.theta)*(obj.strain(1) - obj.strain(2))/(obj.strain(1)+obj.strain(2));
            strainNorm = norm(obj.strain);
            phi = atan2(obj.strain(1)/strainNorm,obj.strain(2)/strainNorm);            
            q3 = ((obj.mu1 + obj.k0)/(obj.k1 - obj.k0) + obj.theta)/tan(phi + pi/4);
            ax = obj.computeAx(q);
            ay = obj.computeAy(q);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = 1 - cParams.volume;
            obj.strain = cParams.strain;
            obj.cx  = cParams.cx;
            obj.cy  = cParams.cy;      
            obj.mu0 = cParams.mu0;
            obj.mu1 = cParams.mu1;
            obj.k0  = cParams.k0;
            obj.k1  = cParams.k1;
        end

        function B = computeB(obj)
            B = obj.mu1*(obj.strain(2) - obj.strain(1) + 2*1i*obj.strain(3));
        end
        
        function C = computeC(obj)
            C = (obj.mu1*(obj.k0 - obj.k1)*(obj.strain(1) + obj.strain(2)))/(obj.mu1+(1-obj.theta)*obj.k0+obj.theta*obj.k1);
        end               
        
        function ax = computeAx(obj,q)
            %ax = (1 + obj.theta + q)/4;
            ax = (1 + obj.theta + q)*obj.cx/2;
        end
        
        function ay = computeAy(obj,q)
            %ay = (1 + obj.theta - q)/4;
            ay = (1 + obj.theta - q)*obj.cy/2;
        end        
        

    end
    
end
