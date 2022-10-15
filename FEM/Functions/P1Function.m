classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
        geometry % !!
    end

    properties (Access = private)
        type
        connec
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nF,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = obj.connec(:,iNode);
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
                end
            end

        end

        function gradFun = computeGradient(obj, quad, mesh)
            % Previous requirements
            obj.createGeometry(quad, mesh);

            % On to calculations
            dNdx = obj.geometry.dNdx;
            nDimf = obj.ndimf;
            nDims = size(dNdx, 1); % derivX, derivY (mesh-related?)
            nNode = size(dNdx, 2);
            nElem = size(dNdx, 3);
            nGaus = size(dNdx, 4);
            
            grad = zeros(nDims,nDimf, nElem, nGaus);
            for iGaus = 1:nGaus
                dNdx_g = dNdx(:,:,:,iGaus);
                for iDims = 1:nDims
                    for iNode = 1:nNode
                        dNdx_i = squeeze(dNdx_g(iDims, iNode,:));
                        nodes = obj.connec(:,iNode);
                        f = obj.fValues(nodes,:);
                        p = (dNdx_i.*f)';
                        pp(1,:,:) = p;
                        grad(iDims,:,:,:) = grad(iDims,:,:,:) + pp;
                    end
                end
            end
            s.fValues    = reshape(grad, [4, 1,6400]);
%             s.ndimf      = nDimf;
            s.quadrature = quad;
            gradFun = FGaussDiscontinuousFunction(s);
        end

        function symGradFun = computeSymmetricGradient(obj, quad, mesh)
            % Previous requirements
            obj.createGeometry(quad, mesh);
            
            % On to calculations
            d_u = obj.convertFValuesToColumn();
            obj.interpolation.computeShapeDeriv(quad.posgp);
            conn = obj.connec;
            nStre   = obj.getNstre();
            nElem   = size(conn,1);
            nNodeEl = size(conn,2);
            nUnkn   = obj.ndimf;
            nGaus   = quad.ngaus;
            symGrad = zeros(nStre,nElem,nGaus);
            for igaus = 1:nGaus
                Bmat = obj.computeB(igaus);
                for istre=1:nStre
                    for inode=1:nNodeEl
                        nodes = conn(:,inode);
                        for idime = 1:nUnkn
                            dofs = nUnkn*(nodes - 1) + idime;
                            idof = nUnkn*(inode - 1) + idime;
                            disp(idof)
                            B = squeeze(Bmat(istre,idof,:));
                            u = d_u(dofs);
                            symGrad(istre,:,igaus) = symGrad(istre,:,igaus) + (B.*u)';
                        end
                        
                    end
                end
            end
            s.fValues    = permute(symGrad, [1 3 2]);
            s.ndimf      = nStre;
            s.quadrature = quad;
            symGradFun = FGaussDiscontinuousFunction(s);
        end

        function plot(obj, m) % 2D domains only
            x = m.coord(:,1);
            y = m.coord(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = obj.fValues(:,idim);
                a = trisurf(m.connec,x,y,z);
                view(0,90)
    %             colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function Bmat = computeB(obj,igaus)
            d.nvoigt = obj.getNstre();
            d.nnodeElem = size(obj.connec,2);
            d.ndofsElem = obj.ndimf*(d.nnodeElem);
            d.ndimf = obj.ndimf;
            s.dim          = d;
            s.geometry     = obj.geometry;
            Bcomp = BMatrixComputer(s);
            Bmat = Bcomp.compute(igaus);
        end

        function createGeometry(obj, quad, mesh)
            int = obj.interpolation;
            int.computeShapeDeriv(quad.posgp);
            s.mesh = mesh;
            g = Geometry.create(s);
            g.computeGeometry(quad,int);
            obj.geometry = g;
        end

        function nstre = getNstre(obj)
            nstreVals = [2, 3, 6];
            nstre = nstreVals(obj.ndimf);
        end

        function fCol = convertFValuesToColumn(obj)
            nVals = size(obj.fValues,1)*size(obj.fValues,2);
            fCol = zeros(nVals,1);
            for idim = 1:obj.ndimf
                fCol(idim:obj.ndimf:nVals, 1) = obj.fValues(:,idim)';
            end
        end

    end

end