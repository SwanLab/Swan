classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        unfittedMesh
    end

    properties (Access = private)
        ndimf
        fun
    end

    methods (Access = public)

        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
        end

        function f = obtainFunctionAtExternalBoundary(obj,iBoundary)
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

            % Provisional solution: ----------------------------
            coord      = obj.unfittedMesh.backgroundMesh.coord;
            diffc      = diff(coord);
            coord(:,all(diffc == 0))=[];
            bCoord     = s.uMesh.backgroundMesh.coord;
            bCoord(:,all(diffc == 0))=[];
            if size(coord,2) == 1
                newFValues = interp1(coord,f,bCoord);
            else
                F          = scatteredInterpolant(coord,f);
                newFValues = F(bCoord);
            end

            ss.fValues = newFValues;
            ss.mesh    = s.uMesh.backgroundMesh;
            s.fun   = P1Function(ss); % ------------------------

             % PENDING - The solution that should be:-----------
             % fNew2        = obj.unfittedMesh.obtainFunctionAtCutMesh(f);


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

            % Provisional solution:------------------------
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
            fNew.fValues = newFValues; % ----------------------

            % PENDING - The solution that should be:-----------
%             fNew2        = obj.unfittedMesh.obtainFunctionAtCutMesh(f);

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

    end
end