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

        function fxV = evaluate(obj, xV)
            assert(isequal(xV, obj.quadrature.posgp), 'Gauss points do not match')
            fxV = obj.fValues;
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
            pp1.mesh   = obj.mesh;
            pp1.connec = obj.mesh.connec;
            projP1 = Projector_toP1(pp1);
            p1fg = projP1.project(obj);
            p1fg.plot();
        end

        function print(obj, s)
            s.mesh = obj.mesh;
            s.fun  = {obj};
            p = FunctionPrinter(s);
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