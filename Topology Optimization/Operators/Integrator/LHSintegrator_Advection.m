classdef LHSintegrator_Advection < LHSintegrator
    
    properties (Access = private)
        geometry
        b
        quadType
        field
    end

    methods (Access = public)
        
        function obj = LHSintegrator_Advection(cParams)
            obj.mesh         = cParams.mesh;
            obj.globalConnec = cParams.globalConnec;
            obj.b            = cParams.b;
            obj.quadType     = cParams.quadType;    
            obj.field        = cParams.field;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function [CX,CY,DX,DY,EX,EY] = compute(obj)
            [CX,CY,DX,DY,EX,EY] = obj.computeElementalLHS();
            CX = obj.assembleMatrixField(CX);
            CY = obj.assembleMatrixField(CY);
            DX = obj.assembleMatrixField(DX);
            DY = obj.assembleMatrixField(DY);
            EX = obj.assembleMatrixField(EX);
            EY = obj.assembleMatrixField(EY);
        end
        
    end
    
   methods (Access = protected)
        
        function [CX,CY,DX,DY,EX,EY] = computeElementalLHS(obj)
            dvolu = obj.mesh.computeDvolume(obj.quadrature);
            ngaus = obj.quadrature.ngaus;
            nelem = obj.mesh.nelem;
    %        nstre = obj.dim.nstre;
            ndpe  = obj.field.dim.ndofsElem;
    %        lhs = zeros(ndpe,ndpe,nelem);
            dN = obj.geometry.dNdx;
            N  = obj.interpolation.shape;

            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = obj.b;
            bF = FeFunction(s);
            bElem = bF.fElem;
            cX = zeros(ndpe,ndpe,nelem);
            cY = zeros(ndpe,ndpe,nelem);
            dX = zeros(ndpe,ndpe,nelem);
            dY = zeros(ndpe,ndpe,nelem);    
            eX = zeros(ndpe,ndpe,nelem);
            eY = zeros(ndpe,ndpe,nelem);  

            for igaus = 1:ngaus
                dNg = dN(:,:,:,igaus);
                dV(:,1) = dvolu(igaus, :);                
                for iNode = 1:3
                    dNi = squeezeParticular(dNg(:,iNode,:),2);
                    Ni = N(iNode,igaus);
                    for jNode = 1:3
                        dNj = squeezeParticular(dNg(:,jNode,:),2);   
                        Nj = N(jNode,igaus);
                        for kNode = 1:3
                            bXk = squeeze(bElem(1,kNode,:));
                            bYk = squeeze(bElem(2,kNode,:));
                            Nk = N(kNode,igaus);
                            dNk = squeezeParticular(dNg(:,kNode,:),2);   
                            

                            dNxi(:,1) = squeeze(dNi(1,:));
                            dNyi(:,1) = squeeze(dNi(2,:));
                            dNxj(:,1) = squeeze(dNj(1,:));
                            dNyj(:,1) = squeeze(dNj(2,:));
                            dNxk(:,1) = squeeze(dNk(1,:));
                            dNyk(:,1) = squeeze(dNk(2,:));

                            
                            ck = Ni.*(dNxj.*dNxk + dNyj.*dNyk);
                            dk = Nk.*(dNxj.*dNxi + dNyj.*dNyi);

                            cXv(1,1,:) = bYk.*ck.*dV;
                            dXv(1,1,:) = bYk.*dk.*dV;

                            cX(iNode,jNode,:) = cX(iNode,jNode,:) + cXv;
                            dX(iNode,jNode,:) = dX(iNode,jNode,:) + dXv;

                            cYv(1,1,:) = bXk.*ck.*dV;
                            dYv(1,1,:) = bXk.*dk.*dV;

                            cY(iNode,jNode,:) = cY(iNode,jNode,:) + cYv;
                            dY(iNode,jNode,:) = dY(iNode,jNode,:) + dYv;
                            

                            e = Ni.*Nj.*Nk;
                            eXV(1,1,:) = e.*bXk.*dV;
                            eYV(1,1,:) = e.*bYk.*dV;

                            eX(iNode,jNode,:) = eX(iNode,jNode,:) + eXV;
                            eY(iNode,jNode,:) = eY(iNode,jNode,:) + eYV;

                        end                   
                    end
                end
            end
            CX = cX;
            CY = cY;
            DX = dX;
            DY = dY;
            EX = eX;
            EY = eY;
        end


       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature(obj.quadType);
           obj.quadrature = quad;
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

        function LHS = assembleMatrixField(obj, lhs)
            s.dim          = obj.field.dim;
            s.globalConnec = obj.field.connec;
            s.nnodeEl      = obj.field.dim.nnodeElem;
            assembler = Assembler(s);
            LHS = assembler.assembleFields(lhs, obj.field, obj.field);
        end        
       
   end
    
end