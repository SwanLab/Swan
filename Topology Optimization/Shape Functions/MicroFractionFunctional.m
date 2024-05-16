classdef MicroFractionFunctional < handle

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
        function obj = MicroFractionFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterField(xD);
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            J  = obj.computeFunction();
            dJ = obj.computeGradient(dC);
            dJ = obj.filterField(dJ);
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

        function xR = filterField(obj,x)
            xR = obj.filter.compute(x,2);
        end

        function J = computeFunction(obj)
            [JAB,JBA,JAA,JBB] = obj.computeAllAlphaBetaValues();
            J = JAB/JAA + JBA/JBB;
        end

        function dJ = computeGradient(obj,dC)
            [JAB,JBA,JAA,JBB] = obj.computeAllAlphaBetaValues();
            a     = obj.alpha;
            b     = obj.beta;
            op    = obj.computeChomogGradientOperation(dC);
            beta1 = (JAA*b - JAB*a)/(JAA^2);
            beta2 = (JBB*a - JBA*b)/(JBB^2);
            dJ1   = obj.computeGradientOfInverse(op,a,beta1);
            dJ2   = obj.computeGradientOfInverse(op,b,beta2);
            dJ    = dJ1+dJ2;
        end

        function [JAB,JBA,JAA,JBB] = computeAllAlphaBetaValues(obj)
            Ch  = obj.stateProblem.Chomog;
            a   = obj.alpha;
            b   = obj.beta;
            JAB = a'*(Ch\b);
            JBA = b'*(Ch\a);
            JAA = a'*(Ch\a);
            JBB = b'*(Ch\b);
        end

        function dChOp = computeChomogGradientOperation(obj,dC)
            tstrain = obj.stateProblem.strainFun;
            dChOp   = @(xV) [];
            for i = 1:length(tstrain)
                opj = @(xV) [];
                for j = 1:length(tstrain)
                    dStrj  = DDP(dC,tstrain{j});
                    dChVij = -DDP(tstrain{i},dStrj);
                    dChEv  = @(xV) reshape(dChVij.evaluate(xV),1,1,size(xV,2),[]);
                    opj    = @(xV) cat(2,opj(xV),dChEv(xV));
                end
                dChOp = @(xV) cat(1,dChOp(xV),opj(xV));
            end
        end

        function dChInv = computeGradientOfInverse(obj,dChOp,a,b)
            Ch          = obj.stateProblem.Chomog;
            wInv        = (Ch\a)*(b'/Ch);
            s.operation = @(xV) squeezeParticular(sum(wInv.*dChOp(xV),[1,2]),1);
            dChInv      = DomainFunction(s);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'ChomogFraction';
        end
    end
end