classdef InternalForceIntegrator < RHSintegrator
    
    methods (Access = public)
        
        function obj = InternalForceIntegrator(cParams)
            obj.init(cParams)
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end
        
        function Fi = compute(obj, stress, testFunc)
            FiElem = obj.computeElementalInternalForce(stress, testFunc);
            Fi     = obj.assembleIntegrand(FiElem, testFunc);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function Fi_el = computeElementalInternalForce(obj, stress, testFunc)
            sig    = stress.evaluate(obj.quadrature.posgp);
            sig    = obj.revertVoigtNotation(sig);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            dNdx   = testFunc.evaluateCartesianDerivatives(obj.quadrature.posgp);
            nDim   = size(dNdx,1);
            nNode  = size(dNdx,2);
            nGaus  = size(dNdx,3);
            nElem  = size(dNdx,4);
            Fi_el  = zeros(nNode*nDim,nElem);
            for iGaus = 1:nGaus
                sig_scaled = sig.*dV(:,iGaus);
                dN         = squeeze(dNdx(:,:,iGaus,:));
                fi_int     = pagemtimes(sig_scaled,dN);
                Fi_el      = Fi_el + reshape(fi_int,[nNode*nDim,nElem]);
            end
        end

    end

    methods (Static, Access = private)

        function Fi = assembleIntegrand(FiElem, test)
            integrand = pagetranspose(FiElem);
            connec    = test.getConnec();
            nDofs     = max(max(connec));
            nDofElem  = size(connec,2);
            Fi        = zeros(nDofs,1);
            for idof = 1:nDofElem
                int = integrand(:,idof);
                con = connec(:,idof);
                Fi  = Fi + accumarray(con,int,[nDofs,1],@sum,0);
            end
        end

        function f_new = revertVoigtNotation(f)
            ndim = size(f,1);
            nel  = size(f,3);
            switch ndim
                case 1
                    f_new = f;
                case 3
                    f_new        = zeros(2,2,nel);
                    f_new(1,1,:) = f(1,:,:);
                    f_new(2,2,:) = f(2,:,:);
                    f_new(1,2,:) = f(3,:,:)/2;
                    f_new(2,1,:) = f_new(1,2,:);
                case 6
                    f_new        = zeros(3,3,nel);
                    f_new(1,1,:) = f(1,:,:);
                    f_new(2,2,:) = f(2,:,:);
                    f_new(3,3,:) = f(3,:,:);
                    f_new(1,2,:) = f(6,:,:)/2;
                    f_new(2,1,:) = f_new(1,2,:);
                    f_new(1,3,:) = f(5,:,:)/2;
                    f_new(3,1,:) = f_new(1,3,:);
                    f_new(2,3,:) = f(4,:,:)/2;
                    f_new(3,2,:) = f_new(2,3,:);
            end
        end

    end
    
end