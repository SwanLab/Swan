classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        fun
    end

    methods (Access = public)

        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
        end

        function f = obtainFunctionAtExternalBoundary(obj,iBoundary) % Only for FeFun so far...
            meshes  = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
            s.uMesh = meshes{iBoundary};
            fType   = obj.fun.fType;


            switch fType
                case 'L2' % Provisional
                    fP1 = obj.fun.project('P1');
                    f   = fP1.fValues;
                otherwise
                    f   = obj.fun.fValues;
            end

            coord      = obj.unfittedMesh.backgroundMesh.coord;
            F          = scatteredInterpolant(coord,f);
            bCoord     = s.uMesh.backgroundMesh.coord;
            ss.fValues = F(bCoord);
            ss.mesh    = s.uMesh.backgroundMesh;


            s.fun   = P1Function(ss);
            f       = UnfittedFunction(s);
        end

        function fxV = evaluateCutElements(obj,q)
            funClass = obj.fun.fType;
            switch funClass
                case 'L2'
                    f = obj.fun.project('P1');
                case 'FE'
                    f = obj.fun;
            end
            c            = obj.unfittedMesh.backgroundMesh.coord;
            c            = c(:,sum(diff(c),1)~=0);
            meshNew      = obj.unfittedMesh.boundaryCutMesh.mesh;
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
            xV = q.posgp;
            fxV = fNew.evaluate(xV);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh   = cParams.uMesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

    end
end