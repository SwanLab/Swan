classdef LHSintegrator_Advection < LHSintegrator
    
    properties (Access = private)
        geometry
        b
    end

    methods (Access = public)
        
        function obj = LHSintegrator_Advection(cParams)
            obj.init(cParams);
            obj.b = cParams.b;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function [LHSX,LHSY] = compute(obj)
            [lhsX,lhsY] = obj.computeElementalLHS();
            LHSX = obj.assembleMatrix(lhsX);
            LHSY = obj.assembleMatrix(lhsY);
        end
        
    end
    
   methods (Access = protected)
        
        function [lhsX,lhsY] = computeElementalLHS(obj)
            dvolu = obj.mesh.computeDvolume(obj.quadrature);
            ngaus = obj.quadrature.ngaus;
            nelem = obj.mesh.nelem;
            nstre = obj.dim.nstre;
            ndpe  = obj.dim.ndofPerElement;
            lhs = zeros(ndpe,ndpe,nelem);
            dN = obj.geometry.dNdx;
            N  = obj.interpolation.shape;

            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = obj.b;
            bF = FeFunction(s);
            bElem = bF.fElem;
            aX = zeros(ndpe,ndpe,nelem);
            aY = zeros(ndpe,ndpe,nelem);

            for igaus = 1:ngaus
                dNg = dN(:,:,:,igaus);
                dV(:,1) = dvolu(igaus, :);                
                for iNode = 1:3
                    dNi = squeezeParticular(dNg(:,iNode,:),2);
                    for jNode = 1:3
                        dNj = squeezeParticular(dNg(:,jNode,:),2);                       
                        for kNode = 1:3
                            bXk = squeeze(bElem(1,kNode,:));
                            bYk = squeeze(bElem(2,kNode,:));
                            Nk = N(kNode,igaus);

                            dNxi(:,1) = squeeze(dNi(1,:));
                            dNyi(:,1) = squeeze(dNi(2,:));
                            dNxj(:,1) = squeeze(dNj(1,:));
                            dNyj(:,1) = squeeze(dNj(2,:));

                            intX(1,1,:) = dNxi.*(dNyj.*Nk.*bXk - dNxj.*Nk.*bYk).*dV;
                            aX(iNode,jNode,:) = aX(iNode,jNode,:) + intX;

                            intY(1,1,:) = dNyi.*(dNyj.*Nk.*bXk - dNxj.*Nk.*bYk).*dV;
                            aY(iNode,jNode,:) = aY(iNode,jNode,:) + intY;                            

                        end                   
                    end
                end
            end
            lhsX = aY;
            lhsY = aY;
        end
        
   end
    
   methods (Access = private)
       
        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function Bcomp = createBComputer(obj)
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = obj.globalConnec;
            s.dofsInElem   = obj.dofsInElem;
            Bcomp = BMatrixComputer(s);
        end
       
   end
    
end