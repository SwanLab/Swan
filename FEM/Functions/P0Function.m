classdef P0Function < FeFunction
    
    properties (Access = public)
    end
    
    properties (Access = private)
        interpolation
    end

    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createFvaluesByElem();
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            % Its a p0 function, so no true need to interpolate -- the
            % value is constant
            nGaus = size(xV,2);
            nFlds = size(obj.fValues,1);
            nElem = size(obj.fValues,3);
            fxV = zeros(nFlds,nGaus,nElem);
            for iGaus = 1:nGaus
                fxV(:,iGaus,:) = squeeze(obj.fValues);
            end
        end

        function dofConnec = computeDofConnectivity(obj)
            conne  = obj.mesh.connec;
            nDimf  = obj.ndimf;
            nDofsE = nDimf;
            nElem  = size(conne,1);
            dofsElem  = zeros(nDofsE,nElem);
            for iUnkn = 1:nDimf
                dofsElem(iUnkn,:) = iUnkn:nDimf:(nElem*nDimf);
            end
            %a = dofsElem(:);
            %dofsElem = reshape(a',[],2)';
            dofConnec = dofsElem;
        end
        
        function N = computeShapeFunctions(obj, quad)
            xV = quad.posgp;
            obj.interpolation.computeShapeDeriv(xV);
            N = obj.interpolation.shape;
        end

        function plot(obj)
            p1DiscFun = obj.project('P1D');
            p1DiscFun.plot();
        end

        function print(obj, s)
            s.mesh = obj.mesh;
            s.fun = {obj};
            p = FunctionPrinter.create(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            nElem = size(obj.mesh.connec, 1);
            nGaus = q.ngaus;

            s.nDimf   = obj.ndimf;
            s.nData   = nElem*nGaus;
            s.nGroup  = nElem;
            s.fValues = obj.getFormattedFValues();
            fps = FunctionPrintingSettings(s);
            [res, pformat] = fps.getDataToPrint();
        end
    end

    methods (Access = public, Static)

        function p0 = create(mesh, ndimf)
            s.fValues = zeros(mesh.nelem,ndimf);
            s.mesh    = mesh;
            p0 = P0Function(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.mesh    = cParams.mesh;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = 'LINEAR';                        
        end

        function createInterpolation(obj)
            m.type = obj.mesh.type;
            obj.interpolation = Interpolation.create(m,'CONSTANT');
        end

        function createFvaluesByElem(obj)
            f = obj.fValues;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fValues = reshape(f',[nDime, 1, nElem]);
        end

        % Printing
        function fM = getFormattedFValues(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            fV = obj.evaluate(q.posgp);
            nGaus   = q.ngaus;
            nComp   = obj.ndimf;
            nElem   = size(obj.mesh.connec, 1);
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