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
            xD = x.obtainDomainFunction();
            xR = obj.filterField(xD);
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            J  = obj.computeFunction();
            dJ = obj.computeGradient(dC{1});
            dJ = obj.filterField({dJ});
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J             = obj.computeNonDimensionalValue(J);
            dJ{1}.fValues = obj.computeNonDimensionalValue(dJ{1}.fValues);
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
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function J = computeFunction(obj)
            Ch = obj.stateProblem.Chomog;
            a  = obj.alpha;
            b  = obj.beta;
            J  = a'*(Ch\b);
        end

        function dJ = computeGradient(obj,dC)
            op = obj.computeChomogGradientOperation(dC);
            dJ = obj.computeGradientOfInverse(op);
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

        function dChInv = computeGradientOfInverse(obj,dChOp)
            Ch          = obj.stateProblem.Chomog;
            a           = obj.alpha;
            b           = obj.beta;
            wInv        = (Ch\a)*(b'/Ch);
            f           = @(xV) squeezeParticular(sum(wInv.*dChOp(xV),[1,2]),1);
            dChInv      = DomainFunction.create(f,obj.mesh);
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