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
            dJ = obj.computeGradient(dC{1});
            dJ = obj.filterField({dJ});
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J     = obj.computeNonDimensionalValue(J);
            dJVal = obj.computeNonDimensionalValue(dJ{1}.fValues);
            dJ{1}.setFValues(dJVal);
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
            Sh  = inv4D(Ch);
            a   = obj.alpha;
            b   = obj.beta;
            JAB = obj.doubleDDP(a,Sh,b);
            JBA = obj.doubleDDP(b,Sh,a);
            JAA = obj.doubleDDP(a,Sh,a);
            JBB = obj.doubleDDP(b,Sh,b);
        end

        function dChOp = computeChomogGradientOperation(obj,dC)
            tstrain = obj.stateProblem.strain;
            dChOp   = @(xV) 0;
            for i = 1:length(tstrain)
                for j = 1:length(tstrain)
                    dStrj  = DDP(dC,tstrain{j});
                    dChVij = -DDP(tstrain{i},dStrj);
                    dChEv  = @(xV) reshape(dChVij.evaluate(xV),1,1,1,1,size(xV,2),[]);
                    eij    = obj.defineCanonicalTensor(i,j);
                    dChOp  = @(xV) dChOp(xV) + dChEv(xV).*eij;
                end
            end
        end

        function v = computeBasesPosition(obj)
            switch obj.mesh.ndim
                case 2
                    v = [1,1; 2,2; 1,2];
                case 3
                    v = [1,1; 2,2; 3,3; 2,3; 1,3; 1,2];
            end
        end

        function eij = defineCanonicalTensor(obj,i,j)
            v   = obj.computeBasesPosition();
            vI  = v(i,:);
            vJ  = v(j,:);
            dim = obj.mesh.ndim;
            eij = zeros(dim,dim,dim,dim);
            eij(vI(1),vI(2),vJ(1),vJ(2)) = 1;
            eij(vI(1),vI(2),vJ(2),vJ(1)) = 1;
            eij(vI(2),vI(1),vJ(1),vJ(2)) = 1;
            eij(vI(2),vI(1),vJ(2),vJ(1)) = 1;
        end

        function dChInv = computeGradientOfInverse(obj,dChOp,a,b)
            Ch     = obj.stateProblem.Chomog;
            Sh     = inv4D(Ch);
            Shb    = tensorprod(Sh,b,[3,4],[1,2]);
            dCShb  = @(xV) tensorprod(dChOp(xV),Shb,[3,4],[1,2]);
            f      = @(xV) reshape(obj.doubleDDP(a,Sh,dCShb(xV)),1,size(xV,2),[]);
            dChInv = DomainFunction.create(f,obj.mesh);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end
    
    methods (Static, Access = private)
        function J = doubleDDP(a,S,b)
            Sb = tensorprod(S,b,[3,4],[1,2]);
            J  = tensorprod(a,Sb,[1,2],[1,2]);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'ChomogFraction';
        end
    end
end