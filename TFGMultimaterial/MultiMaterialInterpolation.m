classdef MultiMaterialInterpolation < handle

   properties (Access = private)
        youngVec
        nuVec
        ndim
   end

   properties (Access = private)
       simpAlls
       muRef
       kappaRef
   end

    methods (Access = public)
        function obj = MultiMaterialInterpolation(cParams)
            obj.init(cParams);
            obj.computeReferenceShearBulk();
            obj.createSimpalls();
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,x)
            nMat       = length(obj.youngVec);
            [mu,kappa] = obj.simpAlls{nMat-2,nMat-1}.computeConsitutiveTensor(x{end});
            for i = nMat-3:-1:1
                m          = obj.createSimpall(obj.muRef{i},obj.kappaRef{i},mu,kappa);
                [mu,kappa] = m.computeConsitutiveTensor(x{i+1});
            end
            mAll       = obj.createSimpall(obj.muRef{end},obj.kappaRef{end},mu,kappa);
            [mu,kappa] = mAll.computeConsitutiveTensor(x{1});
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,x)
            m            = x{1}.mesh;
            I            = LagrangianFunction.create(m,1,'P1');
            I.fValues(:) = 1;
            Z            = LagrangianFunction.create(m,1,'P1');
            [dmu,dkappa] = computeTopologicalDerivatives(obj,I,Z);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.youngVec = cParams.E;
            obj.nuVec    = cParams.nu;
            obj.ndim     = cParams.ndim;
        end

        function computeReferenceShearBulk(obj)
            E  = obj.youngVec;
            nu = obj.nuVec;
            mu = cell(length(E),1);
            k  = cell(length(E),1);
            for i = 1:length(E)
                mu{i} = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E(i),nu(i));
                k{i}  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E(i),nu(i),obj.ndim);
            end
            obj.muRef    = mu;
            obj.kappaRef = k;
        end

        function createSimpalls(obj)
            for i = 1:length(obj.youngVec)-1
                for j = 2:length(obj.youngVec)
                    if i<j
                        obj.simpAlls{i,j} = obj.createSimpall(obj.muRef{i},obj.kappaRef{i},obj.muRef{j},obj.kappaRef{j});
                    end
                end
            end
        end

        function m = createSimpall(obj,mu0,kappa0,mu1,kappa1)
            s.interpolation = 'SIMPALL';
            s.dim           = obj.computeNDimChar();
            s.matA.shear    = mu0;
            s.matA.bulk     = kappa0;
            s.matB.shear    = mu1;
            s.matB.bulk     = kappa1;
            m               = MaterialInterpolator.create(s);
        end

        function dim = computeNDimChar(obj)
            switch obj.ndim
                case 2
                  dim = '2D';
                case 3'
                  dim = '3D';
            end
        end

        function [dmu,dkappa] = computeTopologicalDerivatives(obj,I,Z)
            nMat   = length(obj.youngVec);
            dmu    = cell(nMat,nMat);
            dkappa = cell(nMat,nMat);
            for i = 1:nMat
                for j = 1:nMat
                    if i==j
                        dmu{i,j}    = Z;
                        dkappa{i,j} = Z;
                    elseif i<j
                        [dmu{i,j},dkappa{i,j}] = obj.simpAlls{i,j}.computeConsitutiveTensorDerivative(Z);
                    else
                        [dmu{i,j},dkappa{i,j}] = obj.simpAlls{j,i}.computeConsitutiveTensorDerivative(I);
                        dmu{i,j}               = -dmu{i,j};
                        dkappa{i,j}            = -dkappa{i,j};
                    end
                end
            end
        end
    end
end