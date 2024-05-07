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
            tstrain = obj.stateProblem.strainFun;
            op = @(xV) [];
            for i = 1:length(tstrain)
                opj = @(xV) [];
                for j = 1:length(tstrain)
                    dStrj  = DDP(dC,tstrain{j});
                    dChVij = -DDP(tstrain{i},dStrj);
                    dChEv  = @(xV) reshape(dChVij.evaluate(xV),1,1,size(xV,2),[]);
                    opj    = @(xV) cat(2,opj(xV),dChEv(xV));
                end
                op = @(xV) cat(1,op(xV),opj(xV));
            end
            wInv = Sh*a*b'*Sh;
            s.operation = @(xV) squeezeParticular(sum(wInv.*op(xV),[1,2]),1);
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