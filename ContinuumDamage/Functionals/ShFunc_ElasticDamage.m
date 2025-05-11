classdef ShFunc_ElasticDamage < handle

    properties (Access = private)
        material
        mesh
        quadOrder
    end

    properties (Access = private)
        test
        RHS

        rOld
        r
        r0
        r1

        qLaw
    end

    methods (Access = public)

        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams);
            obj.rOld = copy(obj.r0);
            obj.r    = copy(obj.r0);
            obj.defineRHSIntegrator();
        end

        function setTestFunction(obj,u)
            obj.test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
        end

        function [energy,C] = computeFunction(obj,u)
            obj.computeDamage();
            C = obj.material.obtainTensorSecant();
            e  = SymGrad(u);
            s  = DDP(e,C);
            en = DDP(s,e);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;
        end

        function res = computeResidual(obj,u)
            obj.computeDamage();
            C = obj.material.obtainTensorSecant();
            epsi = SymGrad(u);
            stress = DDP(epsi,C);
            res = obj.RHS.compute(stress,obj.test);
        end

        function [K,Ktan] = computeDerivativeResidual(obj,u,isLoading)
            obj.computeDamage();
            %Ksec = obj.computeDerivativeResidualSecant();
            Ktan = obj.computeDerivativeResidualTangent(u,isLoading);
            K = Ktan;
        end

        function computeDamageEvolutionParam(obj,u)
            C = obj.material.obtainNonDamagedTensor();
            epsi = SymGrad(u);
            tauEpsilon = sqrt(DDP(DDP(epsi,C),epsi));
            tauEpsilon = project(tauEpsilon,obj.r.order);

            fV = zeros(size(obj.r.fValues));

            nodesNoDamage = tauEpsilon.fValues <= obj.rOld.fValues;
            fV(nodesNoDamage) = obj.rOld.fValues(nodesNoDamage);
            fV(~nodesNoDamage) = tauEpsilon.fValues(~nodesNoDamage);

            obj.r.setFValues(fV);

            obj.qLaw.isDamaging = @(xV) tauEpsilon.evaluate(xV) - obj.rOld.evaluate(xV) <= 0 ;
            obj.qLaw.isOverR1 = @(xV) obj.r.evaluate(xV) > obj.r1.evaluate(xV);
        end

        function setROld (obj)
            obj.rOld.setFValues(obj.r.fValues());
        end

        function d = getDamage(obj)
            obj.computeDamage();
            d = obj.material.getDamage();
        end

        function r = getR(obj)
            r = obj.r;
        end

        function qFun = getQ(obj)
            qFun = obj.material.getQ();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.quadOrder = cParams.quadOrder;
            obj.material  = cParams.material;
            obj.r0        = cParams.r0;
            obj.r1        = cParams.r1;

        end

        function LHS = createElasticLHS(obj,material)
            s.type = 'ElasticStiffnessMatrix';
            s.quadratureOrder = obj.quadOrder;
            s.mesh = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.test;
            s.material = material;
            LHS = LHSIntegrator.create(s);
        end
        
        function defineRHSIntegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder= obj.quadOrder;
            s.mesh = obj.mesh;
            obj.RHS = RHSIntegrator.create(s);
        end


        function computeDamage(obj)
            obj.updateQlaw();
            obj.material.updateMaterial(obj.qLaw);
        end

        function sec = computeDerivativeResidualSecant(obj)
            mat = obj.material.obtainTensorSecant();
            LHS = obj.createElasticLHS(mat);
            sec = LHS.compute();
        end

        function deffineQlaw (obj,cParams)
            obj.qLaw.type = cParams.qLaw.type;
            obj.qLaw.A    = cParams.qLaw.A;
            obj.qLaw.H    = cParams.H;
            obj.qLaw.r    = cParams.r;
            obj.qLaw.r0   = cParams.r0;
            obj.qLaw.r1   = cParams.r1;
            obj.qLaw.qInf = cParams.qLaw.qInf;
            obj.qLaw.isDamaging = @(xV) obj.r.evaluate(xV).*0;
            obj.qLaw.isOverR1   = @(xV) obj.r.evaluate(xV).*0;
        end

        function updateQlaw (obj)
            obj.qLaw.r = obj.r;
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,isLoading) 
            if (isLoading)
                C = obj.material.obtainTensorTanget(u); %Ctan
            else
                C = obj.material.obtainTensorSecant(); %Csec
            end
            LHS = obj.createElasticLHS(C);
            tan = LHS.compute();
        end

                % function defineLHSIntegrator(obj)
        %     s.type = 'ElasticStiffnessMatrix';
        %     s.quadratureOrder = obj.quadOrder;
        %     s.mesh  = obj.mesh;
        %     s.test  = obj.test;
        %     s.trial = obj.test;
        %     obj.LHS = LHSIntegrator.create(s);
        % end


        % function rState = getStateR (obj,r,r0,r1)
        %     rState.rHigherR0 = @(xV) (r(xV) > r0(xV));
        %     rState.rHigherR1 = @(xV) (r(xV) > r1(xV));
        % end

        % function q = computeHardening(obj)
        %     switch (obj.qLaw)
        %         case 'Linear'
        %             qInf = @(r, r0, r1) (r0 + obj.H * (r1 - r0));
        %             qLoading =  @(r, r0, r1) (r0 + obj.H .* (r - r0));
        %             q = @(r, r0, r1) (r > r1) .* qInf(r,r0,r1) + (r <= r1) .* qLoading(r,r0,r1);
        %         case 'Exp'
        %             q = @(r, r0, r1) obj.qInfExp -(obj.qInfExp-r0)*exp(obj.A*(1-r/r0));
        % 
        %         otherwise
        %             msg = "Q law is not deffined, put either 'Linear' or 'Exp";
        %             error(msg);
        %     end
        % end
        % 
        % function qDeriv = computeHardeningDerivative(obj)
        %     switch (obj.qLaw)
        %         case 'Linear'
        %             qDeriv = @(r,r0) obj.H;
        %         case 'Exp'
        %             qDeriv = @(r,r0) obj.A*((obj.qInfExp-r0)/r0)*exp(obj.A*(1-r/r0));
        %         otherwise
        %             msg = "Q law is not defined, put either 'Linear' or 'Exp";
        %             error(msg);
        %     end
        % end


        % function tan = computeDerivativeResidualTangent(obj,u,isLoading)
        %     Csec = obj.material.obtainTensor(obj.d);
        %     if (isLoading)
        %         C = obj.material.obtainNonDamagedTensor();
        %         epsi = SymGrad(u);
        %         sigBar = DDP(epsi,C);
        %         q = obj.computeHardening();
        %         qDeriv = obj.computeHardeningDerivative();
        % 
        %         r = @(xV) obj.r.evaluate(xV);
        %         rOld = @(xV) obj.rOld.evaluate(xV);
        %         r0 = @(xV) obj.r0.evaluate(xV);
        %         r1 = @(xV) obj.r1.evaluate(xV);
        % 
        %         rHigherR0 = @(xV) (r(xV) > r0(xV));
        %         rHigherR1 = @(xV) (r(xV) > r1(xV).*isequal(obj.qLaw , 'Linear'));
        % 
        % 
        %         d_dotLoading = @(xV)  (q(r(xV),r0(xV),r1(xV))-qDeriv(r(xV),r0(xV))*r(xV))./(r(xV).^3);
        %         d_dotQInf = @(xV)     (q(r(xV),r0(xV),r1(xV)))./(r(xV).^3);
        % 
        %         checkDamageLimit = @(xV) rHigherR1(xV).*d_dotQInf(xV);
        %         checkLoading = @(xV) d_dotLoading(xV).*(rHigherR0(xV).*~rHigherR1(xV));
        % 
        %         op = @(xV)  checkDamageLimit(xV) + checkLoading(xV);
        % 
        %         d_dot = DomainFunction.create(op,obj.mesh);
        %         Ctan2 = Expand(d_dot).*OP(sigBar,sigBar);
        % 
        %         op = @(xV) Csec.evaluate(xV) - Ctan2.evaluate(xV);
        %         Ctan = DomainFunction.create(op,obj.mesh);
        %     else
        %         Ctan = Csec;
        %     end
        %     LHS = obj.createElasticLHS(Ctan);
        %     tan = LHS.compute();
        % end

    end
end