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
            Ch  = obj.stateProblem.Chomog;
            a   = obj.alpha;
            b   = obj.beta;
            JAB = a'*(Ch\b);
            JBA = b'*(Ch\a);
            JAA = a'*(Ch\a);
            JBB = b'*(Ch\b);
            J   = JAB/JAA + JBA/JBB;
        end

        function dJ = computeGradient(obj,dC)
            Ch      = obj.stateProblem.Chomog;
            a       = obj.alpha;
            b       = obj.beta;
            JAB     = a'*(Ch\b);
            JBA     = b'*(Ch\a);
            JAA     = a'*(Ch\a);
            JBB     = b'*(Ch\b);
            tstrain = obj.stateProblem.strainFun;
            op      = @(xV) [];
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
            beta1       = (JAA*b - JAB*a)/(JAA^2);
            beta2       = (JBB*a - JBA*b)/(JBB^2);
            wInv1       = (Ch\a)*(beta1'/Ch);
            wInv2       = (Ch\b)*(beta2'/Ch);
            s.operation = @(xV) squeezeParticular(sum(wInv1.*op(xV),[1,2]),1);
            dJ1         = DomainFunction(s);
            s.operation = @(xV) squeezeParticular(sum(wInv2.*op(xV),[1,2]),1);
            dJ2         = DomainFunction(s);
            dJ          = dJ1+dJ2;
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