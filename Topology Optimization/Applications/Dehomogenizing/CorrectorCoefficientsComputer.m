classdef CorrectorCoefficientsComputer < handle

    properties (Access = private)
        quadrature
        dVolum
        oCorrectorDerivative
        LHS
        RHS
    end

    properties (Access = private)
        mesh
        orthogonalCorrector
    end

    methods (Access = public)

        function obj = CorrectorCoefficientsComputer(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createDvolum();
        end

        function cOpt = compute(obj,b)
            obj.createOrthogonalCorrectorDerivatives();
            obj.computeLHS();
            obj.computeRHS(b);
            c = obj.computeReferenceCoefficients();
            cOpt = obj.computeBestIntegerCoeff(c);
            % cOpt = c;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                = cParams.mesh;
            obj.orthogonalCorrector = cParams.orthogonalCorrector;
        end

        function cOpt = computeBestIntegerCoeff(obj,c)
            Jref = obj.computeCost(c);
            w    = obj.computeAllBinaryPossibilities(c);
            Jall = obj.computeCostOfAllIntegerSimilarCoeff(w,c);
            [Jopt,iOpt] = min(Jall);
            wOpt = w(iOpt,:);
            cOpt = obj.computeIntegerCoeff(wOpt,c);
            %cOpt = c;
        end

        function w = computeAllBinaryPossibilities(obj,c)
            nSing = length(c);
            w  = ff2n(nSing);
        end

        function J = computeCostOfAllIntegerSimilarCoeff(obj,w,c)
            nComb  = size(w,1);
            J    = zeros(nComb,1);
            for iComb = 1:nComb
                wI = w(iComb,:);
                cI = obj.computeIntegerCoeff(wI,c);
                J(iComb) = obj.computeCost(cI);
            end
        end

        function cI = computeIntegerCoeff(obj,wI,c)
            cFloor = floor(c);
            cCeil  = ceil(c);
            cI = (1-wI').*cFloor + wI'.*cCeil;
        end

        function J = computeCost(obj,c)
            A = obj.LHS;
            b = obj.RHS;
            r = (A*c - b);
            J = 0.5*(r')*r;
        end

        function c = computeReferenceCoefficients(obj)
            c = obj.LHS\obj.RHS;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            obj.quadrature = q;
        end

        function createOrthogonalCorrectorDerivatives(obj)
            psi = obj.orthogonalCorrector;
            xV  = obj.quadrature.posgp;
            nSing = numel(psi);
            dP = cell(nSing,1);
            for iSing = 1:nSing
                dPsiV = psi{iSing}.evaluateGradient(xV);
                dP{iSing} = dPsiV.fValues;
            end
            obj.oCorrectorDerivative = dP;
        end

        function createDvolum(obj)
            q = obj.quadrature;
            nDim = obj.mesh.ndim;
            dV(1,:,:) = obj.mesh.computeDvolume(q);
            dV = repmat(dV,nDim,1,1);
            obj.dVolum = dV;
        end

        function computeLHS(obj)
            dP  = obj.oCorrectorDerivative;
            dV  = obj.dVolum;
            nSing = numel(obj.orthogonalCorrector);
            lhs = zeros(nSing,nSing);
            for iS = 1:nSing
                dPi = dP{iS};
                for jS = 1:nSing
                    dPj = dP{jS};
                    lhsIJ   = dPi.*dPj.*dV;
                    lhs(iS,jS) = sum(lhsIJ(:));
                end
            end
            obj.LHS = lhs;
        end

        function computeRHS(obj,b)
            bG    = obj.computeOrientationInGauss(b);
            nSing = numel(obj.orthogonalCorrector);
            rhs = zeros(nSing,1);
            dP = obj.oCorrectorDerivative;
            dV   = obj.dVolum;
            for iS = 1:nSing
                dPi = dP{iS};
                rhsI = dPi.*bG.*dV;
                rhs(iS) = sum(rhsI(:));
            end
            obj.RHS = rhs;
        end

        function bfG = computeOrientationInGauss(obj,b)
            q      = obj.quadrature;
            xGauss = q.posgp;
            bfG    = b.evaluate(xGauss);
            %bfG    = permute(b.evaluate(xGauss),[1 3 2]);
        end

    end

end