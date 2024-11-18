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
            I   = LagrangianFunction.create(rho{1}.mesh,1,'P1');
            I.fValues(:) = 1;
            Z   = LagrangianFunction.create(rho{1}.mesh,1,'P1');

            [dmu14,dkappa14] = obj.simpAlls.m14.computeConsitutiveTensorDerivative(I);
            dkSeba = dkappa{4,1};
            error1 = (dkSeba - dkappa14) ./ dkappa14;
            error1 = error1.project('P1',rho{1}.mesh);
            error1 = sum(error1.fValues);

            dmSeba = dmu{4,1};
            error2 = (dmSeba - dmu14) ./ dmu14;
            error2 = error2.project('P1',rho{1}.mesh);
            error2 = sum(error2.fValues);


            [dmu12,dkappa12] = obj.simpAlls.m12.computeConsitutiveTensorDerivative(Z);
            ratio2 = (   (10/3 * obj.youngVec(2)*dkappa{2,1})   ./    dkappa12   );
            ratio2 = ratio2.project('P1',rho{1}.mesh);
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
            nu    = 0.25; % !!!
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
            chiVal         = obj.splitCellIntoValues(chi);
            C              = obj.computeC(chiVal);
            Cm             = obj.elasticTensor{j};
            
            P              = obj.computeGradientCoefficientsMatrix(chi,i,j);
            dC             = pagemtimes(C,P);
            dCij           = dC;
        end

        function Cv = computeC(obj,chiVal)
            C      = obj.elasticTensor{1};            
            nElem  = size(chiVal,2);            
            tgamma = obj.computeTgamma3(chiVal);
            Cv     = repmat(C,[1 1 nElem]);
            Cv     = Cv.*tgamma;                        
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