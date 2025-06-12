classdef LocalDamageFunctional < handle
    
    properties (Access = private)
        mesh
        dissipation
        %constant
        constantFun
        l0
        testPhi
    end
    
    methods (Access = public)
        
        function obj = LocalDamageFunctional(cParams)
            obj.init(cParams)
        end
        
        function F = computeFunctional(obj,phi,quadOrder)
            alphaFun = obj.obtainDissipationFunction(phi,'Function');

            int = Integrator.create('Function',obj.mesh,quadOrder);
            %F   = (obj.constant/obj.l0)*int.compute(alphaFun);
            constant = obj.updateConstant(phi);
            F   = (1/obj.l0)*int.compute(alphaFun.*constant);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            dAlphaFun =  obj.obtainDissipationFunction(phi,'Jacobian');

            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSIntegrator.create(s);
            %J   = (obj.constant/obj.l0)*RHS.compute(dAlphaFun,obj.testPhi);
            constant = obj.updateConstant(phi);
            J   = (1/obj.l0)*RHS.compute(dAlphaFun.*constant,obj.testPhi);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            ddAlphaFun =  obj.obtainDissipationFunction(phi,'Hessian');  
            s.trial    = obj.testPhi;
            s.test     = obj.testPhi;
            %s.function = ddAlphaFun;
            s.mesh     = obj.mesh;
            s.type     = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            %LHS = LHSIntegrator.create(s);
            %H = (obj.constant/obj.l0)*LHS.compute();
            constant = obj.updateConstant(phi);
            s.function = ddAlphaFun.*constant;
            LHS = LHSIntegrator.create(s);
            H = (1/obj.l0)*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.testPhi     = copy(cParams.testSpace.phi);          
            obj.dissipation = cParams.dissipation.interpolation;
            %obj.constant    = cParams.dissipation.constant;
            obj.l0          = cParams.l0;
        end

        function c = updateConstant(obj,phi)
            isOverThreshold = (phi.fun>=0.25);
            c = (isOverThreshold).*(0.5*(0.9^2)*0.1) + (~isOverThreshold).*(0.5*(1^2)*0.1);
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