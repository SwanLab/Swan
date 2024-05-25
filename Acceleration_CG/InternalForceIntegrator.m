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
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            dNdx   = testFunc.evaluateCartesianDerivatives(obj.quadrature.posgp);

            nDim   = size(dNdx,1);
            nNode  = size(dNdx,2);
            nGaus  = size(dNdx,3);
            nElem  = size(dNdx,4);
            Fi_el  = zeros(nNode*nDim,nElem);

            dNdx  = obj.transformToVoigt(dNdx);
            for iGaus = 1:nGaus
                sig_scaled = sig.*dV(:,iGaus);
                fi_int     = pagemtimes(squeeze(dNdx(:,:,iGaus,:)),sig_scaled);
                Fi_el      = Fi_el + squeeze(fi_int(:,1,:));
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

        function f_new = transformToVoigt(f)
            ndim   = size(f,1);
            nShape = size(f,2);
            nGauss = size(f,3);
            nel    = size(f,4);
            switch ndim
                case 1
                    f_new = f;
                case 2 % xx - yy - xy
                    f_new = zeros(nShape*ndim,3,nGauss,nel);
                    for iShape = 1:nShape
                        f_new(2*iShape-1,[1 3],:,:) = [f(1,iShape,:,:), f(2,iShape,:,:)];
                        f_new(2*iShape,[2,3],:,:)   = [f(2,iShape,:,:), f(1,iShape,:,:)];
                    end
                case 3 % xx - yy - xy - xz - yz
                    f_new = zeros(nShape*ndim,6,nGauss,nel);
                    for iShape = 1:nShape
                        f_new(3*iShape-2,[1 4 5],:,:) = [f(1,iShape,:,:), f(2,iShape,:,:), f(3,iShape,:,:)];
                        f_new(3*iShape-1,[2 4 6],:,:) = [f(2,iShape,:,:), f(1,iShape,:,:), f(3,iShape,:,:)];
                        f_new(3*iShape,[3 5 6],:,:)   = [f(3,iShape,:,:), f(1,iShape,:,:), f(2,iShape,:,:)];
                        
                    end
            end
        end

    end
    
end