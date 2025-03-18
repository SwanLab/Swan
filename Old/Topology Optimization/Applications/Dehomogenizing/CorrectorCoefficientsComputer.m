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
        end

        function cOpt = compute(obj,b)
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

        function computeLHS(obj)
            psi  = obj.orthogonalCorrector;
            nSing = numel(obj.orthogonalCorrector);
            lhs = zeros(nSing,nSing);
            for iS = 1:nSing
                psiI = psi{iS};                                
                for jS = 1:nSing
                    psiJ = psi{jS};
                    lhs(iS,jS) = Integrate(DDP(Grad(psiI),Grad(psiJ)));
                end
            end
            obj.LHS = lhs;
        end

        function computeRHS(obj,b)
            nSing = numel(obj.orthogonalCorrector);
            rhs = zeros(nSing,1);
            for iS = 1:nSing
                psiS    = obj.orthogonalCorrector{iS};                
                rhs(iS) = Integrate(DDP(Grad(psiS),b));
            end
            obj.RHS = rhs;
        end

    end

end