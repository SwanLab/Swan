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
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,x,C)
            for i = 1:length(obj.youngVec)
                for j = 1:length(obj.youngVec)
                    dC      = obj.computeTensorDerivativeIJ(x{1},C,i,j);
                    dmuVal  = squeeze(dC(3,3,:));
                    dlamVal = squeeze(dC(1,2,:));
                    dmu{i,j}     = obj.computeP1Function(x{1}{1}.mesh,dmuVal);
                    s.operation = @(xV) dmu{i,j}.evaluate(xV);
                    dmu{i,j}    = DomainFunction(s);
                    dlam    = obj.computeP1Function(x{1}{1}.mesh,dlamVal);
                    N       = x{1}{1}.mesh.ndim;
                    dkappa{i,j}  = obj.computeBulkMagnitude(dlam,dmu{i,j},N);
                end
            end





            rho = x{2};
            m = rho{1}.mesh;
            I   = LagrangianFunction.create(m,1,'P1');
            I.fValues(:) = 1;
            Z   = LagrangianFunction.create(m,1,'P1');

            [dmu12,dkappa12] = obj.simpAlls.m12.computeConsitutiveTensorDerivative(I);
            [dmu13,dkappa13] = obj.simpAlls.m13.computeConsitutiveTensorDerivative(I);
            [dmu14,dkappa14] = obj.simpAlls.m14.computeConsitutiveTensorDerivative(I);
            [dmu23,dkappa23] = obj.simpAlls.m23.computeConsitutiveTensorDerivative(I);
            [dmu24,dkappa24] = obj.simpAlls.m24.computeConsitutiveTensorDerivative(I);
            [dmu34,dkappa34] = obj.simpAlls.m34.computeConsitutiveTensorDerivative(I);

            ratio12K = Mean(    dkappa{1,2}./dkappa12   ,m,2);
            ratio13K = Mean(    dkappa{1,3}./dkappa13   ,m,2);
            ratio14K = Mean(    dkappa{1,4}./dkappa14   ,m,2);
            ratio23K = Mean(    dkappa{2,3}./dkappa23   ,m,2);
            ratio24K = Mean(    dkappa{2,4}./dkappa24   ,m,2);
            ratio34K = Mean(    dkappa{3,4}./dkappa34   ,m,2);

            ratio12m = Mean(    dmu{1,2}./dmu12         ,m,2);
            ratio13m = Mean(    dmu{1,3}./dmu13         ,m,2);
            ratio14m = Mean(    dmu{1,4}./dmu14         ,m,2);
            ratio23m = Mean(    dmu{2,3}./dmu23         ,m,2);
            ratio24m = Mean(    dmu{2,4}./dmu24         ,m,2);
            ratio34m = Mean(    dmu{3,4}./dmu34         ,m,2); % All fine but different sign


            [dmu21,dkappa21] = obj.simpAlls.m12.computeConsitutiveTensorDerivative(Z);
            ratio21K = Mean(    dkappa{2,1}./dkappa21   ,m,2);
            ratio21m = Mean(    dmu{2,1}./dmu21         ,m,2); % Here all fine, even the sign
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

        function coefMatrix2 = computeGradientCoefficientsMatrix(obj,chi,i,j)
            chiVal = obj.splitCellIntoValues(chi);
            nu    = 1/3; % !!!
            alpha  = (1+nu)/(1-nu);
            beta = (3-nu)/(1+nu);
            E     = obj.youngVec;
            E1    = E*chiVal;
            %i     = obj.currentMat;
            %j     = obj.inclusionMat;
            g     = E(j)/E(i);
            c1    = -(1+beta)/(1+g*beta)*(1-g)./ones(size(E1));
            c2    = -0.5*(alpha-beta)/(1+g*beta)*(g*(g-2)+1)/(1+g*alpha)./ones(size(E1));

            coefMatrix2(1,1,:) = c1+c2;
            coefMatrix2(1,2,:) = c2;
            coefMatrix2(1,3,:) = 0;
            coefMatrix2(2,1,:) = c2;
            coefMatrix2(2,2,:) = c1+c2;
            coefMatrix2(2,3,:) = 0;
            coefMatrix2(3,1,:) = 0;
            coefMatrix2(3,2,:) = 0;
            coefMatrix2(3,3,:) = c1;
        end

        function dCij = computeTensorDerivativeIJ(obj,chi,C2,i,j)
            Cm     = obj.elasticTensor{i};
            P      = obj.computeGradientCoefficientsMatrix(chi,i,j);
            dC     = pagemtimes(Cm,P);
            dCij   = dC;
        end

        function tgamma3 = computeTgamma3(obj,chiVal)
            tgamma         = (obj.youngVec/obj.youngVec(1))*chiVal;
            tgamma3(:,1,:) = [tgamma;tgamma;tgamma];
            tgamma3        = repmat(tgamma3,[1,3,1]);
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

            obj.simpAlls.m12 = obj.createSimpall(E2,nu2,E1,nu1);
            obj.simpAlls.m13 = obj.createSimpall(E3,nu3,E1,nu1);
            obj.simpAlls.m14 = obj.createSimpall(E4,nu4,E1,nu1);
            obj.simpAlls.m23 = obj.createSimpall(E3,nu3,E2,nu2);
            obj.simpAlls.m24 = obj.createSimpall(E4,nu4,E2,nu2);
            obj.simpAlls.m34 = obj.createSimpall(E4,nu4,E3,nu3);
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
    end
end