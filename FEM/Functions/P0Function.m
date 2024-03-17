classdef P0Function < FeFunction
    
    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
    end
    
    properties (Access = private)
        interpolation
    end

    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.computeNDofs();
        end

        function fxV = evaluate(obj, xV)
            % Its a p0 function, so no true need to interpolate -- the
            % value is constant
            nGaus = size(xV,2);
            nFlds = size(obj.fValues,2);
            nElem = size(obj.fValues,1);
            fxV = zeros(nFlds,nGaus,nElem);
            fVals = reshape(obj.fValues',[nFlds, 1, nElem]);
            for iGaus = 1:nGaus
                fxV(:,iGaus,:) = squeeze(fVals);
            end
        end

        function [p0sub, mesh_sub, l2g] = restrictTo(obj, domain)
            [mesh_sub, l2g] = obj.mesh.getBoundarySubmesh(domain);
            dofs = sum(ismember(obj.mesh.connec, l2g(mesh_sub.connec)), 2) == 2;
            s.fValues = obj.fValues(dofs, :);
            s.mesh    = mesh_sub;
            p0sub = P0Function(s);
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
        
        function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
        end

        function plot(obj)
            p1Fun = obj.project('P1');
            p1Fun.plot();
            % p1Fun.plotLine();
        end

        function plotArrowVector(obj)
            %figure()
            a = obj.fValues;
            xy = obj.mesh.computeBaricenter();
            x = xy(1,:)';
            y = xy(2,:)';
            ax = squeeze(a(1,:,:));
            ay = squeeze(a(2,:,:));
            q = quiver(x,y,ax,ay);
            q.ShowArrowHead = 'off';
        end        

        function print(obj, filename, software)
            if nargin == 2; software = 'GiD'; end
            s.mesh = obj.mesh;
            s.fun = {obj};
            s.type = software;
            s.filename = filename;
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
            obj.order   = 'CONSTANT';
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,'CONSTANT');
        end

        function computeNDofs(obj)
            obj.nDofsElem = obj.ndimf*obj.interpolation.nnode;
            obj.nDofs = obj.ndimf * size(obj.fValues, 1);
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
