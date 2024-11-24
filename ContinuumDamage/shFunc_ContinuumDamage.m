classdef shFunc_ContinuumDamage < handle

    properties (Access = public)

    end
    properties (Access = private)
        material
        mesh

        internalDamage
        internalElastic
        external

        r0
        H
    end

    methods (Access = public)

        function obj = shFunc_ContinuumDamage(cParams)
            obj.init(cParams);
        end

        function totalEnergy = computeTotalEnergyDamage (obj, quadOrder, u,r,fext)
            internalEnergy = obj.computeInternalFunction(quadOrder,u,r);
            externalEnergy = obj.computeExternalFunction(u,fext,quadOrder);
            totalEnergy = internalEnergy - externalEnergy;
        end

        % EXTERNAL WORK FUNCTIONS
        function F = computeExternalFunction(obj,u,fExt,quadOrder)
            bMesh = obj.mesh.createBoundaryMesh{4}; %CARA SUPERIOR
            int = Integrator.create('Function',bMesh.mesh,quadOrder);
            [u, fExt] = obj.adaptFuns(u,fExt);
            F = int.compute(u.*fExt);
        end

        function Ju = computeGradient(obj,u,fExt,quadOrder)
            bMesh = obj.mesh.createBoundaryMesh{4}; %CARA SUPERIOR
            s.mesh = bMesh.mesh;
            s.quadType = quadOrder;
            s.type = 'ShapeFunction';
            RHS = RHSintegrator.create(s);

            [u,fExt] = adaptFuns(obj,u,fExt);
            test = LagrangianFunction.create(bMesh.mesh,u.ndimf,u.order);
            Ju = RHS.compute(fExt,test);
            Ju = obj.reducedToFull(Ju,bMesh);
        end
        % INTERNAL ENERGY FUNCTIONS

        function energy = computeInternalFunction(obj,quadOrder,u,r)

            d = obj.computeDamage(r);
            Cdamage = obj.material.obtainTensor(d);

            epsi = SymGrad(u);
            funct = DDP(DDP(epsi,Cdamage),epsi);
            energy = 0.5*(Integrator.compute(funct,obj.mesh,quadOrder));

        end

        function jacobian = computeJacobian(obj,quadOrder,u,r)

            d = obj.computeDamage(r);
            Cdamage = obj.material.obtainTensor(d);
            epsi = SymGrad(u);
            b = DDP(epsi,Cdamage);

            S.type = 'ShapeSymmetricDerivative';
            S.quadratureOrder=quadOrder;
            S.mesh = obj.mesh;

            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            rhs = RHSintegrator.create(S);
            jacobian = rhs.compute(b,test);

        end

        function hessian = computeHessian(obj,quadOrder,u,r)

            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            d = obj.computeDamage(r);

            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = 2;
            S.mesh = obj.mesh;
            S.material = obj.material.obtainTensor(d);
            S.test = test;
            S.trial = test;

            lhs = LHSintegrator.create(S);
            hessian = lhs.compute();
        end

        function rOut = newState (obj,rIn,u)

            C = obj.material.obtainNonDamagedTensor;
            epsi = SymGrad(u);
            tauEpsi = power(DDP(DDP(epsi,C),epsi),0.5);

            if tauEpsi <= rIn

                rOut = rIn;
            else
                rOut = tauEpsi;
            end
        end

        function d = computeDamage(obj,r)
            q = obj.computeHardening();
            s.operation = @(xV) 1-(q.evaluate(r.evaluate(xV))./(r.evaluate(xV)));
            s.ndimf = 1;
            d = DomainFunction(s);
        end
    end

    methods (Access = private)
        function init (obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r0 = cParams.r0;
            obj.H = cParams.H;
        end

        function [uFun,fExtFun] = adaptFuns(obj,u,fExt) %% ADAPTING TO TOP BOUNDARY %%
            bMesh = obj.mesh.createBoundaryMesh{4};
            nodes = unique(bMesh.globalConnec);
            if isempty(fExt)
                uFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                uFun.fValues = u.fValues(nodes,:);
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
            else
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                fExtFun.fValues = fExt.fValues(nodes,:);
                uFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                uFun.fValues = u.fValues(nodes,:);
            end
        end

        function JuFull = reducedToFull(obj,Ju,bMesh)
            nNodes = obj.mesh.nnodes;
            nDim = obj.mesh.ndim;
            nDofs = nNodes*obj.mesh.ndim;
            JuFull = zeros(nDofs,1);
            nodes = unique(bMesh.globalConnec);
            nNodesB = length(nodes);
            ForceDofs = zeros(nNodesB*nDim,1);
            for iDim = 1:nDim
                ForceDofs(iDim:nDim:(end-nDim+iDim)) = nDim*(nodes-1)+iDim;
            end
            JuFull(ForceDofs) = Ju;
        end

        function q = computeHardening(obj)

            s.operation = @(r) obj.r0 + obj.H *(r-obj.r0);
            s.ndimf = 1;
            q = DomainFunction(s);

        end


    end
end