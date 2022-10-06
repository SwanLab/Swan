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

        function strain = computeGradientStrain(obj, quad, mesh)
            % Previous requirements
            obj.createGeometry(quad, mesh);
            
            % On to calculations
            d_u(1:2:29,1) = obj.fValues(:,1)';
            d_u(2:2:30,1) = obj.fValues(:,2)';

            obj.interpolation.computeShapeDeriv(quad.posgp);
            conn = obj.connec;
            nStre   = obj.getNstre();
            nElem   = size(conn,1);
            nNodeEl = size(conn,2);
            nUnkn   = obj.ndimf;
            nGaus   = quad.ngaus;
            strain = zeros(nStre,nElem,nGaus);
            for igaus = 1:nGaus
                Bmat = obj.computeB(igaus);
                for istre=1:nStre
                    for inode=1:nNodeEl
                        nodes = conn(:,inode);
                        for idime = 1:nUnkn
                            dofs = nUnkn*(nodes - 1) + idime;
                            ievab = nUnkn*(inode-1)+idime;
                            B = squeeze(Bmat(istre,ievab,:));
                            u = d_u(dofs);
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(B.*u)';
                        end
                        
                    end
                end
            end
            strain = permute(strain, [3 1 2]);
        end

        function plot(obj, m) % 2D domains only
            dim = 1;
            x = m.coord(:,1);
            y = m.coord(:,2);
            z = obj.fValues(:,dim);
            figure()
            a = trisurf(m.connec,x,y,z);
            view(0,90)
%             colorbar
            shading interp
            a.EdgeColor = [0 0 0];
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
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
            Bmat = Bcomp.computeBmat(igaus);
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

    end

end