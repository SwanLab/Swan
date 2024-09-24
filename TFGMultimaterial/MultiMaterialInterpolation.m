classdef MultiMaterialInterpolation < handle

   properties (Access = private)
        multiYoung
        C1
        nu1
   end

    methods (Access = public)
        function obj = MultiMaterialInterpolation(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,x)
            chiFun     = x.obtainDomainFunction;
            [~,chi]    = chiFun.computeAtNodesAndElements();
            tgamma     = (obj.multiYoung/obj.multiYoung(1))*chi.fValues';
            Ceff       = obj.C1*tgamma;
            lambdaVals = Ceff(6,:);
            muVals     = Ceff(4,:);

            s.order     = 'P0';
            s.fValues   = lambdaVals';
            s.mesh      = x.levelSets{1}.fun.mesh;
            lambda = LagrangianFunction(s);

            s.order   = 'P0';
            s.fValues = muVals';
            s.mesh    = chi.mesh;
            mu        = LagrangianFunction(s);
            
            N     = chi.mesh.ndim;
            kappa = lambda + 2.*mu/N;
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,x)
            chiFun     = x.obtainDomainFunction;
            [~,chi]    = chiFun.computeAtNodesAndElements();
            tgamma     = (obj.multiYoung/obj.multiYoung(1))*chi.fValues';

            nu    = obj.nu1;
            alpha = (3-nu)/(1+nu);
            beta  = (1+nu)/(1-nu);
            E2    = obj.multiYoung(2);
            la0   = nu.*E2./((1+nu).*(1-2.*nu));
            mu0   = E2./(2.*(1+nu));
            la0   = 2.*mu0.*la0./(la0+2.*mu0); % plane stress
            C     = [la0+2*mu0 0 la0; 0 2*mu0 0; la0 0 la0+2*mu0];
            nMat  = length(obj.multiYoung);
            E     = obj.multiYoung;
            tE    = E(1)*tgamma; 
            nElem = length(tE);
            for i = 1:nMat
                for j = 1:nMat
                    g  = E(j)/E(i);
                    c1 = 0.5*((1-g)./(1+alpha*g))./tE;
                    c2 = c1.*((g.*(alpha-2*beta)-1)./(1+beta*g));


                    coefMatrix2(1,1,:) = 4*c1+c2;
                    coefMatrix2(1,2,:) = 0;
                    coefMatrix2(1,3,:) = c2;
                    coefMatrix2(2,1,:) = 0;
                    coefMatrix2(2,2,:) = 8*c1;
                    coefMatrix2(2,3,:) = 0;
                    coefMatrix2(3,1,:) = c2;
                    coefMatrix2(3,2,:) = 0;
                    coefMatrix2(3,3,:) = 4*c1+c2;

                    Cv = repmat(C,[1 1 nElem]);
                    CvDC = pagemtimes(Cv,coefMatrix2);
                    CvDC2 = pagemtimes(CvDC,C);
                    dC(:,:,:,i,j) = CvDC2;

                    dmuVal  = squeeze(dC(3,3,:,i,j));
                    dlamVal = squeeze(dC(1,2,:,i,j));

                    s.order     = 'P0';
                    s.fValues   = dlamVal;
                    s.mesh      = x.levelSets{1}.fun.mesh;
                    dlam = LagrangianFunction(s);

                    s.order   = 'P0';
                    s.fValues = dmuVal;
                    s.mesh    = chi.mesh;
                    dmu{i,j}        = LagrangianFunction(s);
                    N     = chi.mesh.ndim;
                    dkappa{i,j} = dlam + 2.*dmu{i,j}/N;
                end
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.multiYoung = cParams.E;
            obj.C1         = cParams.C0;
            obj.nu1        = cParams.nu1;
        end
    end
end