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

        function fxV = evaluateCutElements(obj,xV)
            % FeFunction:
            %             mesh          = obj.unfittedMesh.backgroundMesh;
            %             coordOriginal = mesh.coord;
            %             n0            = size(coordOriginal,1)+1;
            %             bCMesh        = obj.unfittedMesh.boundaryCutMesh.mesh;
%             coordComplete = bCMesh.coord;
%             oldfValues = obj.fun.fValues;
%             x          = coordOriginal(:,1);
%             y          = coordOriginal(:,2);
%             F          = scatteredInterpolant(x,y,oldfValues);
%             newfValues = F(coordComplete(n0:end,:));
%             s.fValues = [oldfValues;newfValues];
%             s.mesh    = bCMesh;
%             newFun    = eval([class(obj.fun),'(s)']);
%             fxV       = newFun.evaluate(xV);

            % L2:
            %             obj.fun.updateMesh(bCMesh.mesh);
            %             fxV       = obj.fun.evaluate(xV);
            %             obj.fun.updateMesh(mesh);


            % Old:
            fType = obj.fun.fType;
            switch fType
                case 'L2' % Provisional
                    fP1        = obj.fun.project('P1');
                    s.fValues  = fP1.fValues;
                    mesh       = obj.unfittedMesh.backgroundMesh;
                    bcMesh     = obj.unfittedMesh.boundaryCutMesh;
                    connec     = mesh.connec;
                    bcConnec   = connec(bcMesh.cellContainingSubcell,:);
                    mmm.connec = bcConnec;
                    mmm.type   = mesh.type;
                    s.mesh     = mmm;
                    f          = P1Function(s);
                    fxV        = f.evaluate(xV);
                otherwise
                    s.fValues  = obj.fun.fValues;
                    mesh       = obj.unfittedMesh.backgroundMesh;
                    bcMesh     = obj.unfittedMesh.boundaryCutMesh;
                    connec     = mesh.connec;
                    bcConnec   = connec(bcMesh.cellContainingSubcell,:);
                    mmm.connec = bcConnec;
                    mmm.type   = mesh.type;
                    s.mesh     = mmm;
                    f          = P1Function(s);
                    fxV        = f.evaluate(xV);
            end
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