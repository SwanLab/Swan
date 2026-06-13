classdef RigidBodyFunction < BaseFunction

    properties (Access = public)
%                  ndimf
        nbasis
        basisFunctions
    end

    properties (Access = private)
        fvalues
        refPoint
        %         mesh

        fun

    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = RigidBodyFunction(cParams)
            obj.init(cParams)
            obj.computeBasisFunctions();
        end

%         function fxV = evaluateNew(obj, xGLoc)
%             phiU = obj.basisFunctions{1}.evaluate(xGLoc);
%             phiV = obj.basisFunctions{2}.evaluate(xGLoc);
%             phiT = obj.basisFunctions{3}.evaluate(xGLoc);
%             u     = obj.fvalues(1);
%             v     = obj.fvalues(2);
%             theta = obj.fvalues(3);
%             fxV = u*phiU + v*phiV + theta*phiT;
%         end

        function bE = evaluateBasisFunctions(obj,xGloc)
            for i=1:obj.nbasis
                bE{i} = obj.basisFunctions{i}.evaluate(xGloc);
            end
        end

        function R = projectBasisFunctions(obj,interp)
            for i = 1:obj.nbasis
                R{i} = obj.basisFunctions{i}.project(interp);
            end
        end

        function plot(obj)
            p1DiscFun = obj.project('P1D');
            p1DiscFun.plot();
        end

         function RB = restrictBasisToBoundaryMesh(obj,bMesh)
%              nodes      = unique(bMesh.globalConnec(:));
             mesh       = bMesh.mesh;
             RB = RigidBodyFunction.create(mesh,obj.refPoint);
%              s.fValues  = zeros(obj.nbasis,1);
%              s.ndimf    = obj.ndimf;  
%              s.refPoint = obj.refPoint;
% %              functionType = obj.functionType;
%              for i=1:obj.nbasis
%                  basis{i} = obj.basisFunctions{i}.fValues(nodes,:);                 
%              end
%              RB = RigidBodyFunction.create(bMesh.mesh,basis,obj.functionType);
         end
         
    end

    methods (Access = public, Static)

        function RB = create(mesh,refPoint)
            s.ndimf    = mesh.ndim;
            s.fvalues  = zeros(mesh.nnodes, s.ndimf);
            s.mesh     = mesh;
            s.refPoint = refPoint;
            RB = RigidBodyFunction(s);
        end

        function fS = times(f1,f2)
            fS = f1.fValues.*f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = RigidBodyFunction(s);
        end

        function fS = sum(f1,f2)
            fS = f1.fValues+f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = RigidBodyFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fvalues  = cParams.fvalues;
            obj.refPoint = cParams.refPoint;
            obj.mesh     = cParams.mesh;
            obj.ndimf    = cParams.ndimf;
        end

        function computeBasisFunctions(obj)
            switch obj.mesh.ndim
                case 2
                    obj.basisFunctions{1} = obj.computeXTranslationBase2D();
                    obj.basisFunctions{2} = obj.computeYTranslationBase2D();
                    obj.basisFunctions{3} = obj.computeRotationBase2D();
                case 3
                    obj.basisFunctions{1} = obj.computeXTranslationBase3D();
                    obj.basisFunctions{2} = obj.computeYTranslationBase3D();
                    obj.basisFunctions{3} = obj.computeZTranslationBase3D();
                    obj.basisFunctions{4} = obj.computeXRotationBase3D();
                    obj.basisFunctions{5} = obj.computeYRotationBase3D();
                    obj.basisFunctions{6} = obj.computeZRotationBase3D();
            end
            obj.nbasis = length(obj.basisFunctions);
        end

        function f = computeXTranslationBase2D(obj)
            f = ConstantFunction.create([1;0],obj.mesh);
        end

        function f = computeYTranslationBase2D(obj)
            f = ConstantFunction.create([0;1],obj.mesh);
        end

        function f = computeRotationBase2D(obj)
            x0 = obj.refPoint(1);
            y0 = obj.refPoint(2);
            s.fHandle = @(x) [-(x(2,:,:)-y0);x(1,:,:)-x0];
            s.ndimf   = obj.ndimf;
            s.mesh    = obj.mesh;
            f= AnalyticalFunction(s);
        end

        function f = computeXTranslationBase3D(obj)
            f = ConstantFunction.create([1;0;0],obj.mesh);
        end

        function f = computeYTranslationBase3D(obj)
            f = ConstantFunction.create([0;1;0],obj.mesh);
        end

        function f = computeZTranslationBase3D(obj)
            f = ConstantFunction.create([0;0;1],obj.mesh);
        end

        function f = computeXRotationBase3D(obj)
            y0 = obj.refPoint(2);
            z0 = obj.refPoint(3);
            s.fHandle = @(x) [zeros(size(x(1,:,:)));-(x(3,:,:)-z0);x(2,:,:)-y0];
            s.ndimf   = obj.ndimf;
            s.mesh    = obj.mesh;
            f= AnalyticalFunction(s);
        end

        function f = computeYRotationBase3D(obj)
            x0 = obj.refPoint(1);
            z0 = obj.refPoint(3);
            s.fHandle = @(x) [x(3,:,:)-z0;zeros(size(x(1,:,:)));-(x(1,:,:)-x0)];
            s.ndimf   = obj.ndimf;
            s.mesh    = obj.mesh;
            f= AnalyticalFunction(s);
        end

        function f = computeZRotationBase3D(obj)
            x0 = obj.refPoint(1);
            y0 = obj.refPoint(2);
            s.fHandle = @(x) [-(x(2,:,:)-y0);x(1,:,:)-x0;zeros(size(x(1,:,:)))];
            s.ndimf   = obj.ndimf;
            s.mesh    = obj.mesh;
            f= AnalyticalFunction(s);
        end

    end

    methods(Access=protected)

        function fxV = evaluateNew(obj, xGLoc)
            switch obj.mesh.ndim
                case 2
                    phiU = obj.basisFunctions{1}.evaluate(xGLoc);
                    phiV = obj.basisFunctions{2}.evaluate(xGLoc);
                    phiT = obj.basisFunctions{3}.evaluate(xGLoc);
                    u     = obj.fvalues(1);
                    v     = obj.fvalues(2);
                    theta = obj.fvalues(3);
                    fxV = u*phiU + v*phiV + theta*phiT;
                case 3
            end
        end

    end

end