classdef UnfittedFunction < L2Function

    properties (Access = public)
        ndimf
        unfittedMesh
    end

    properties (Access = private)
        fun
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.evaluateInnerElements(xV);
        end

        function fxV = evaluateCutElements(obj,q)
            funClass = obj.fun.fType;
            switch funClass
                case 'L2'
                    f = obj.fun.project('P1');
                case 'FE'
                    f = obj.fun;
            end

            % Provisional solution:------------------------
            c            = obj.unfittedMesh.backgroundMesh.coord;
            diffc        = diff(c);
            c(:,all(diffc == 0))=[];
            meshNew      = obj.unfittedMesh.innerCutMesh.mesh;
            cNew         = meshNew.coord;
            cNew(:,all(diffc == 0))=[];
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
            fNew.fValues = newFValues; % -------------------

            % PENDING - The solution that should be:-----------
            fNew2        = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(f);

            xV  = q.posgp;
            fxV = fNew.evaluate(xV);
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