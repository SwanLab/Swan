classdef LocalDamageFunctional < handle
    
    properties (Access = private)
        mesh
        dissipation
        cw
        Gc
        l0
        testPhi
    end
    
    methods (Access = public)
        
        function obj = LocalDamageFunctional(cParams)
            obj.init(cParams)
        end
        
        function F = computeCost(obj,phi,quadOrder)
            alphaFun = obj.obtainDissipationFunction(phi,'Function');
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = (obj.Gc/(obj.cw*obj.l0))*int.compute(alphaFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            dAlphaFun =  obj.obtainDissipationFunction(phi,'Jacobian');
            J = IntegrateRHS(@(v) (obj.Gc/(obj.cw*obj.l0)).*DP(v,dAlphaFun),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            ddAlphaFun =  obj.obtainDissipationFunction(phi,'Hessian');         
            H = IntegrateLHS(@(u,v) (obj.Gc/(obj.cw*obj.l0)).*ddAlphaFun.*DP(v,u),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.testPhi     = copy(cParams.testSpace.phi);          
            obj.dissipation = cParams.dissipation.interpolation;
            obj.cw          = cParams.dissipation.constant;
            obj.Gc          = cParams.Gc;
            obj.l0          = cParams.l0;
        end

        function disFun = obtainDissipationFunction(obj,phi,type)
            switch type
                case 'Function'
                    fun = obj.dissipation.fun;
                case'Jacobian'
                     fun = obj.dissipation.dfun;
                case 'Hessian'
                    fun = obj.dissipation.ddfun;
            end
            s.operation = @(xV) fun.evaluate(phi.evaluate(xV));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            disFun = DomainFunction(s);
        end

    end
end