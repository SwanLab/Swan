 classdef FGaussDiscontinuousFunction < handle
    % nDimf * nGaus * nElem
    properties (Constant, Access = public)
        fType = 'GAUSSPOINTS'
    end

    properties (Access = public)
        ndimf
        fValues
        quadrature
    end

    properties (Access = private)
        mesh
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = FGaussDiscontinuousFunction(cParams)
            obj.init(cParams)
        end

        function fun = project(obj,target)
            s.mesh          = obj.mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end
        
        function fxV = evaluate(obj, xV)
            % assert(isequal(xV, obj.quadrature.posgp), 'Gauss points do not match')
            fxV = obj.fValues;
        end
        
        function dNdx  = computeCartesianDerivatives(obj, quad)
            assert(isequal(quad,obj.quadrature), 'Quadrature does not match');
            nElem = size(obj.mesh.connec,1);
            nNode = obj.mesh.interpolation.nnode;
            nDime = obj.mesh.interpolation.ndime;
            nGaus = quad.ngaus;
            invJ  = obj.mesh.computeInverseJacobian(quad,obj.mesh.interpolation);
            dShapeDx  = zeros(nDime,nNode,nElem,nGaus);
            for igaus = 1:nGaus
                dShapes = obj.mesh.interpolation.deriv(:,:,igaus);
                for jDime = 1:nDime
                    invJ_JI   = invJ(:,jDime,:,igaus);
                    dShape_KJ = dShapes(jDime,:);
                    dSDx_KI   = bsxfun(@times, invJ_JI,dShape_KJ);
                    dShapeDx(:,:,:,igaus) = dShapeDx(:,:,:,igaus) + dSDx_KI;
                end
            end
            dNdx = dShapeDx;
        end

        function applyVoigtNotation(obj)
            switch obj.ndimf
                case 4
                    obj.applyVoigt2D()
                case 9
                    obj.applyVoigt3D()
            end
        end

        function plot(obj)
            s.mesh = obj.mesh;
            s.projectorType = 'P1D';
            proj = Projector.create(s);
            p1fun = proj.project(obj);
            p1fun.plot();
        end

        function print(obj, filename, software)
            if nargin == 2; software = 'Paraview'; end
            s.mesh = obj.mesh;
            s.fun = {obj};
            s.type = software;
            s.filename = filename;
            p = FunctionPrinter.create(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            nElem = size(obj.fValues, 3);
            nGaus = obj.quadrature.ngaus;
            s.nDimf   = obj.ndimf;
            s.nData   = nElem*nGaus;
            s.nGroup  = nElem;
            s.fValues = obj.getFormattedFValues();
            fps = FunctionPrintingSettings(s);
            [res, pformat] = fps.getDataToPrint();
        end

        function dofConnec = computeDofConnectivity(obj)
            % This assumes that FGaussDiscFun comes from a P1Fun...
            conne  = obj.mesh.connec;
            nDimf  = obj.ndimf;
            nNode  = size(conne, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(conne,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = conne(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

        function v = computeL2norm(obj)
            v = Norm.computeL2(obj.mesh,obj,obj.quadrature);
        end        

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues    = cParams.fValues;
            obj.quadrature = cParams.quadrature;
            obj.ndimf      = size(cParams.fValues,1);
            obj.mesh       = cParams.mesh;
        end

        function applyVoigt2D(obj)
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.fValues,3);
            fV(1,:,:) = obj.fValues(1,:,:); % xx
            fV(2,:,:) = obj.fValues(4,:,:); % yy
            fV(3,:,:) = obj.fValues(2,:,:) + obj.fValues(3,:,:); % xy
            fV = reshape(fV, [3 nGaus nElem]);
            obj.fValues = fV;
            obj.ndimf = 3;
        end

        function applyVoigt3D(obj)
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.fValues,3);
            fV(1,:,:) = obj.fValues(1,:,:); % xx
            fV(2,:,:) = obj.fValues(5,:,:); % yy
            fV(3,:,:) = obj.fValues(9,:,:); % zz
            fV(4,:,:) = obj.fValues(2,:,:) + obj.fValues(4,:,:); % xy
            fV(5,:,:) = obj.fValues(3,:,:) + obj.fValues(7,:,:); % xz
            fV(6,:,:) = obj.fValues(6,:,:) + obj.fValues(8,:,:); % yz
            fV = reshape(fV, [6 nGaus nElem]);
            obj.fValues = fV;
            obj.ndimf = 6;
        end

        % Printing
        function fM = getFormattedFValues(obj)
            fV = obj.fValues;
            nGaus = obj.quadrature.ngaus;
            nComp = obj.ndimf;
            nElem = size(obj.fValues, 3);
            fM  = zeros(nGaus*nElem,nComp);
            for iStre = 1:nComp
                for iGaus = 1:nGaus
                    rows = linspace(iGaus,(nElem - 1)*nGaus + iGaus,nElem);
                    fM(rows,iStre) = fV(iStre,iGaus,:);
                end
            end
        end
        
    end
    
end
