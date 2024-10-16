classdef MultiMaterialInterpolation < handle

   properties (Access = private)
        youngVec
        elasticTensorA
        elasticTensorB
   end

   properties (Access = private)
       currentMat
       inclusionMat
   end

    methods (Access = public)
        function obj = MultiMaterialInterpolation(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,x)
            chi                 = obj.computeElementWiseCharacteristicFunction(x);
            [muVals,lambdaVals] = obj.computeShearLambdaValues(chi);
            mu                  = obj.computeP0Function(chi.mesh,muVals');
            lambda              = obj.computeP0Function(chi.mesh,lambdaVals');            
            N                   = chi.mesh.ndim;
            kappa               = lambda + 2.*mu/N;
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,x)
            chi     = obj.computeElementWiseCharacteristicFunction(x);
            dC      = obj.computeTensorDerivativeIJ(chi);
            dmuVal  = squeeze(dC(2,2,:))/4;
            dlamVal = squeeze(dC(1,3,:));
            dmu     = obj.computeP0Function(chi.mesh,dmuVal);
            dlam    = obj.computeP0Function(chi.mesh,dlamVal);
            N       = chi.mesh.ndim;
            dkappa  = dlam + 2.*dmu/N;
        end

        function computeFromTo(obj,i,j)
            obj.currentMat   = i;
            obj.inclusionMat = j;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.youngVec       = cParams.E;
            obj.elasticTensorA = cParams.CA;
            obj.elasticTensorB = cParams.CB;
        end

        function [muVals,lambdaVals] = computeShearLambdaValues(obj,chi)
            chiVal     = chi.fValues';
            tgamma     = (obj.youngVec/obj.youngVec(1))*chiVal;
            Ceff       = obj.elasticTensorA*tgamma;
            lambdaVals = Ceff(6,:);
            muVals     = Ceff(4,:);
        end

        function coefMatrix2 = computeGradientCoefficientsMatrix(obj,chi)
            nu    = 0.25; % !!!
            alpha = (3-nu)/(1+nu);
            beta  = (1+nu)/(1-nu);
            E     = obj.youngVec;
            E1    = E*chi.fValues';
            i     = obj.currentMat;
            j     = obj.inclusionMat;
            g     = E(j)/E(i);
            c1    = -0.5*((1-g)./(1+alpha*g))./E1;
            c2    = c1.*((g.*(alpha-2*beta)-1)./(1+beta*g));

            coefMatrix2(1,1,:) = 4*c1+c2;
            coefMatrix2(1,2,:) = 0;
            coefMatrix2(1,3,:) = c2;
            coefMatrix2(2,1,:) = 0;
            coefMatrix2(2,2,:) = 8*c1;
            coefMatrix2(2,3,:) = 0;
            coefMatrix2(3,1,:) = c2;
            coefMatrix2(3,2,:) = 0;
            coefMatrix2(3,3,:) = 4*c1+c2;
        end

        function dCij = computeTensorDerivativeIJ(obj,chi)
            nElem          = size(chi.fValues,1);
            C              = obj.elasticTensorB;
            coefMatrix2    = obj.computeGradientCoefficientsMatrix(chi);
            tgamma         = (obj.youngVec/obj.youngVec(1))*chi.fValues';
            tgamma3(:,1,:) = [tgamma;tgamma;tgamma];
            tgamma3        = repmat(tgamma3,[1,3,1]);
            Cv             = repmat(C,[1 1 nElem]);
            CvDC           = pagemtimes(Cv,coefMatrix2);
            dC             = pagemtimes(CvDC,C);
            dCij           = dC.*(tgamma3.^2);
        end
    end

    methods (Static, Access = private)
        function chi = computeElementWiseCharacteristicFunction(x)
            chiFun  = x.obtainDomainFunction;
            [~,chi] = chiFun.computeAtNodesAndElements();
        end

        function f = computeP0Function(m,fValues)
            s.order   = 'P0';
            s.fValues = fValues;
            s.mesh    = m;
            f         = LagrangianFunction(s);
        end
    end
end