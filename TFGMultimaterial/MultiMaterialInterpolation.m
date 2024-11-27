classdef MultiMaterialInterpolation < handle

   properties (Access = private)
        youngVec
        elasticTensor

        simpAlls
   end

    methods (Access = public)
        function obj = MultiMaterialInterpolation(cParams)
            obj.init(cParams);

            obj.createSimpalls(cParams);
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,x)

            [muVals,lambdaVals] = obj.computeShearLambdaValues(x{1});
            mu                  = obj.computeP1Function(x{1}{1}.mesh,muVals);
            lambda              = obj.computeP1Function(x{1}{1}.mesh,lambdaVals);
            N                   = x{1}{1}.mesh.ndim;
            kappa               = obj.computeBulkMagnitude(lambda,mu,N);






            mu0    = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.youngVec(4),1/3);
            kappa0 = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.youngVec(4),1/3,2);

            mu1     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.youngVec(1),1/3);
            kappa1  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.youngVec(1),1/3,2);

            rho1 = x{2}{1};
            rho2 = x{2}{2};
            rho3 = x{2}{3};

            [mu23,kappa23]   = obj.simpAlls{2,3}.computeConsitutiveTensor(rho3);
            m123             = obj.createSimpallFromShearBulk(mu1,kappa1,mu23,kappa23);
            [mu123,kappa123] = m123.computeConsitutiveTensor(rho2);
            mAll             = obj.createSimpallFromShearBulk(mu0,kappa0,mu123,kappa123);
            [mu2,kappa2]       = mAll.computeConsitutiveTensor(rho1);
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,x,C)
            nMat         = length(obj.youngVec);
            m            = x{2}{1}.mesh;
            I            = LagrangianFunction.create(m,1,'P1');
            I.fValues(:) = 1;
            Z            = LagrangianFunction.create(m,1,'P1');
            dmu          = cell(nMat,nMat);
            dkappa       = cell(nMat,nMat);
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

    methods (Access = private)
        function init(obj,cParams)
            obj.youngVec         = cParams.E;
            obj.elasticTensor{1} = cParams.CA;
            obj.elasticTensor{2} = cParams.CB;
            obj.elasticTensor{3} = cParams.CC;
            obj.elasticTensor{4} = cParams.CD;
        end

        function [muVals,lambdaVals] = computeShearLambdaValues(obj,chi)
            chiVal            = obj.splitCellIntoValues(chi);
            tgamma(1,1,:)     = (obj.youngVec/obj.youngVec(1))*chiVal;
            Ceff       = obj.elasticTensor{1}.*tgamma;
            lambdaVals = squeeze(Ceff(1,2,:));
            muVals     = squeeze(Ceff(3,3,:));
        end       

    end

    methods (Static, Access = private)
        function chiVal = splitCellIntoValues(chi)
            chiVal     = zeros(length(chi),length(chi{1}.fValues));
            for i = 1:length(chi)
                chiVal(i,:) = chi{i}.fValues;
            end
        end

        function f = computeP1Function(m,fValues)
            s.order   = 'P1';
            s.fValues = fValues;
            s.mesh    = m;
            f         = LagrangianFunction(s);
        end

        function kappa = computeBulkMagnitude(lambda,mu,N)
            s.operation = @(xV) lambda.evaluate(xV) + 2*mu.evaluate(xV)/N;
            kappa       = DomainFunction(s);
        end
    end






    methods (Access = private)
        function createSimpalls(obj,cParams) % 12, 13, 14, 23, 24, 34
            E1  = cParams.E(1);
            E2  = cParams.E(2);
            E3  = cParams.E(3);
            E4  = cParams.E(4);
            nu1 = cParams.nu(1);
            nu2 = cParams.nu(2);
            nu3 = cParams.nu(3);
            nu4 = cParams.nu(4);

            obj.simpAlls{1,2} = obj.createSimpall(E1,nu1,E2,nu2);
            obj.simpAlls{1,3} = obj.createSimpall(E1,nu1,E3,nu3);
            obj.simpAlls{1,4} = obj.createSimpall(E1,nu1,E4,nu4);
            obj.simpAlls{2,3} = obj.createSimpall(E2,nu2,E3,nu3);
            obj.simpAlls{2,4} = obj.createSimpall(E2,nu2,E4,nu4);
            obj.simpAlls{3,4} = obj.createSimpall(E3,nu3,E4,nu4);
        end

        function m = createSimpall(obj,E0,nu0,E1,nu1)
            ndim = 2;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
        end

        function m = createSimpallFromShearBulk(obj,mu0,kappa0,mu1,kappa1)
            s.interpolation = 'SIMPALL';
            s.dim           = '2D';
            s.matA.shear    = mu0;
            s.matA.bulk     = kappa0;
            s.matB.shear    = mu1;
            s.matB.bulk     = kappa1;
            m               = MaterialInterpolator.create(s);
        end
    end
end