classdef MicroAlphaBetaFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        material
        stateProblem
        alpha
        beta
    end

    methods (Access = public)
        function obj = MicroAlphaBetaFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            obj.material.setDesignVariable(xR);
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            Ch = obj.stateProblem.Chomog;
            Sh = inv(Ch);
            a  = obj.alpha;
            b  = obj.beta;
            J  = a'*Sh*b;
            tstrain = obj.stateProblem.strainFluctFun{1}; % NO HAURIA DE SER STRAINFUN?? eij+epsij
            dStr    = DDP(dC,tstrain);
            dj      = -DDP(tstrain,dStr);
            wInv = Sh*a*b'*Sh;
            s.operation = @(xV) pagemtimes(wInv,dj.evaluate(xV));
            dJ = DomainFunction(s);
            dJ     = obj.filter.compute(dJ,2);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J          = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.material     = cParams.material;
            obj.stateProblem = cParams.stateProblem;
            obj.alpha        = cParams.alpha;
            obj.beta         = cParams.beta;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'ChomogAlphaBeta';
        end
    end
end