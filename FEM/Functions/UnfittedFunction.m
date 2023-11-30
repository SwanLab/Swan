classdef UnfittedFunction < L2Function

    properties (Access = public)
        ndimf
        unfittedMesh
        fun
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.evaluateInnerElements(xV);
        end

        function fxV = evaluateCutElements(obj,xVloc)
            funClass = obj.fun.fType;
            switch funClass
                case 'L2'
                    f = obj.fun.project('P1');
                case 'FE'
                    f = obj.fun;
            end
            c            = obj.unfittedMesh.backgroundMesh.coord;
            c            = c(:,sum(diff(c),1)~=0);
            meshNew      = obj.unfittedMesh.innerCutMesh.mesh;
            cNew         = meshNew.coord;
            cNew         = cNew(:,sum(diff(cNew),1)~=0);
            if size(c,2) == 1
                newFValues = interp1(c,f.fValues,cNew);
            else
                F          = scatteredInterpolant(c,f.fValues);
                newFValues = F(cNew);
            end
            s.feFunType  = class(f);
            s.mesh       = meshNew;
            s.ndimf      = obj.ndimf;
            fNew         = FeFunction.createEmpty(s);
            fNew.fValues = newFValues;
            q = Quadrature.set(meshNew.type);
            q.computeQuadrature('QUADRATIC');
            xV2 = q.posgp;
            fxV = fNew.evaluate(xV2);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh   = cParams.uMesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

        function fxV = evaluateInnerElements(obj,xV)
            fxV     = obj.fun.evaluate(xV);
            gMesh   = obj.unfittedMesh.backgroundMesh;
            inMesh  = obj.unfittedMesh.innerMesh;
            gCoor   = gMesh.computeXgauss(xV);
            inCoor  = inMesh.mesh.computeXgauss(xV);
            gCoor1  = squeeze(gCoor(:,1,:))';
            inCoor1 = squeeze(inCoor(:,1,:))';
            isVoid  = not(ismember(gCoor1,inCoor1,'rows'));
            fxV(:,:,isVoid) = 0;
        end

    end
end