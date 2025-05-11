classdef damageLaw < handle
        
    properties (Access = private)
        qLawClass
        qLawParams

        r
        q
        qDot

        mesh
    end
    
    methods (Access = public)   
        function obj = damageLaw(cParams,mesh)
            obj.init(cParams,mesh)
        end
        
        function updateParams(obj,newValue)

            obj.qLawParams.r = newValue.r;
            obj.qLawParams.isDamaging = newValue.isDamaging;
            obj.qLawParams.isDamaging = newValue.isDamaging;
            obj.qLawParams.isOverR1 = newValue.isOverR1;

            obj.qLawClass = HardeningLaw.create(obj.qLawParams);
            obj.setParams();
        end
        function qFun = getQFun (obj)
            obj.computeQ();
            op = @(xV) obj.q(xV);
            qFun = DomainFunction.create(op,obj.mesh);
        end
        function d = computeDamage(obj)
            op = @(xV) 1-obj.q(xV)/obj.r(xV);
            
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            s.operation = op;

            d = DomainFunction(s);
        end
        
        function dDot = computeDamageDerivative (obj)
            dDot = @(xV) (obj.q(xV) - obj.qDot(xV)*obj.r(xV))/(obj.r(xV)).^3;
        end
    end    
    methods (Access = private)
        function init(obj,cParams,mesh)
            obj.qLawParams = cParams;
            obj.qLawClass = HardeningLaw.create(cParams);
            obj.setParams();
            obj.mesh = mesh;
        end 

        function setParams(obj)
           obj.computeQ();
           obj.computeQDerivative();
           obj.setR();

        end

        function computeQ (obj)
           obj.q = @(xV) obj.qLawClass.computeQ(xV);
        end

        function computeQDerivative(obj)
            obj.qDot = @(xV) obj.qLawClass.computeQDerivative();
        end
        
        function setR(obj)
            obj.r = @(xV) obj.qLawParams.r.evaluate(xV);
        end

        function qDot = computeQderiv (obj)
            qDot = @(xV) obj.qLawClass.computeDerivative();
        end
    end
end