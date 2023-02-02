classdef LHSintegrator_Advection < LHSintegrator
    
    properties (Access = private)
        geometry
        b
        rho
        quadType
        field
    end

    methods (Access = public)
        
        function obj = LHSintegrator_Advection(cParams)
            obj.mesh         = cParams.mesh;
            obj.globalConnec = cParams.globalConnec;
            obj.b            = cParams.b;
            obj.rho          = cParams.rho;
            obj.quadType     = cParams.quadType;    
            obj.field        = cParams.field;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = compute(obj)
            [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = obj.computeElementalLHS();
            CX = obj.assembleMatrixField(CX);
            CY = obj.assembleMatrixField(CY);
            DX = obj.assembleMatrixField(DX);
            DY = obj.assembleMatrixField(DY);
            EX = obj.assembleMatrixField(EX);
            EY = obj.assembleMatrixField(EY);
            Kxx = obj.assembleMatrixField(Kxx);
            Kxy = obj.assembleMatrixField(Kxy);
            Kyy = obj.assembleMatrixField(Kyy);
            Mxx = obj.assembleMatrixField(Mxx);
            Mxy = obj.assembleMatrixField(Mxy);
            Myy = obj.assembleMatrixField(Myy);
            Mrho = obj.assembleMatrixField(Mrho);
        end
        
    end
    
   methods (Access = protected)
        
        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = computeElementalLHS(obj)
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

            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = obj.rho;
            rhoF = FeFunction(s);
            rhoElem = rhoF.fElem;



            cX = zeros(ndpe,ndpe,nelem);
            cY = zeros(ndpe,ndpe,nelem);
            dX = zeros(ndpe,ndpe,nelem);
            dY = zeros(ndpe,ndpe,nelem);    
            eX = zeros(ndpe,ndpe,nelem);
            eY = zeros(ndpe,ndpe,nelem);  
            Kxx = zeros(ndpe,ndpe,nelem);  
            Kyy = zeros(ndpe,ndpe,nelem);  
            Kxy = zeros(ndpe,ndpe,nelem);  
            Mxx = zeros(ndpe,ndpe,nelem);  
            Myy = zeros(ndpe,ndpe,nelem);  
            Mxy = zeros(ndpe,ndpe,nelem);  
            Mrho = zeros(ndpe,ndpe,nelem);  

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
                            rhok = squeeze(rhoElem(1,kNode,:));
                            
                            Nk = N(kNode,igaus);
                            dNk = squeezeParticular(dNg(:,kNode,:),2);   
                            

                            dNidx(:,1) = squeeze(dNi(1,:));
                            dNidy(:,1) = squeeze(dNi(2,:));
                            dNjdx(:,1) = squeeze(dNj(1,:));
                            dNjdy(:,1) = squeeze(dNj(2,:));
                            dNkdx(:,1) = squeeze(dNk(1,:));
                            dNkdy(:,1) = squeeze(dNk(2,:));

                         % version 1   
                            ck = Ni.*(dNjdx.*dNkdx + dNjdy.*dNkdy);
                            dk = Nk.*(dNjdx.*dNidx + dNjdy.*dNidy);

                            cXv(1,1,:) = bYk.*ck.*dV;
                            dXv(1,1,:) = bYk.*dk.*dV;

                            cX(iNode,jNode,:) = cX(iNode,jNode,:) + cXv;
                            dX(iNode,jNode,:) = dX(iNode,jNode,:) + dXv;

                            cYv(1,1,:) = bXk.*ck.*dV;
                            dYv(1,1,:) = bXk.*dk.*dV;

                            cY(iNode,jNode,:) = cY(iNode,jNode,:) + cYv;
                            dY(iNode,jNode,:) = dY(iNode,jNode,:) + dYv;



                            %%% RegularizationTerms
                    
                            for lNode = 1:3
                              

                                Nl = N(lNode,igaus);
                                dNl = squeezeParticular(dNg(:,lNode,:),2);
                                dNldx(:,1) = squeeze(dNl(1,:));
                                dNldy(:,1) = squeeze(dNl(2,:));

                                rhol = squeeze(rhoElem(1,lNode,:));
                                mr(1,1,:) = 4*Ni.*Nj.*(rhok.*Nk).*(1-rhol).*Nl.*dV;
                                Mrho(iNode,jNode,:) = Mrho(iNode,jNode,:) + mr;
                                
                             
                                bXl = squeeze(bElem(1,lNode,:));
                                bYl = squeeze(bElem(2,lNode,:));

                                k0 = (dNidx.*dNjdx + dNidy.*dNjdy).*Nk.*Nl.*dV;
                                kxx(1,1,:) = (bYk.*bYl).*k0;
                                kxy(1,1,:) = (bXk.*bYl).*k0;
                                kyy(1,1,:) = (bXk.*bXl).*k0;

                                m0 = (Ni.*Nj).*(dNldx.*dNkdx+dNldy.*dNkdy).*dV;
                                mxx(1,1,:) = (bYl.*bYk).*m0;
                                mxy(1,1,:) = (bXl.*bYk).*m0;
                                myy(1,1,:) = (bXl.*bXk).*m0;



                             Kxx(iNode,jNode,:) = Kxx(iNode,jNode,:) + kxx;
                             Kxy(iNode,jNode,:) = Kxy(iNode,jNode,:) + kxy;
                             Kyy(iNode,jNode,:) = Kyy(iNode,jNode,:) + kyy;

                             Mxx(iNode,jNode,:) = Mxx(iNode,jNode,:) + mxx;
                             Mxy(iNode,jNode,:) = Mxy(iNode,jNode,:) + mxy;
                             Myy(iNode,jNode,:) = Myy(iNode,jNode,:) + myy;
                            end
 
                         % version 2
% 
%                             cxk = Ni.*dNkdy.*(dNjdx.*bXk + dNjdx.*bYk);
%                             cyk = Ni.*dNkdx.*(dNjdx.*bXk + dNjdx.*bYk);
% 
%                             dxk = Nk.*dNjdx.*(-dNidy.*bXk  + dNidx.*bYk);
%                             dyk = Nk.*dNjdy.*( dNidy.*bXk  - dNidx.*bYk);
% 
%                             cXv(1,1,:) = cxk.*dV;
%                             dXv(1,1,:) = dxk.*dV;
% 
%                             cX(iNode,jNode,:) = cX(iNode,jNode,:) + cXv;
%                             dX(iNode,jNode,:) = dX(iNode,jNode,:) + dXv;
% 
%                             cYv(1,1,:) = cyk.*dV;
%                             dYv(1,1,:) = dyk.*dV;
% 
%                             cY(iNode,jNode,:) = cY(iNode,jNode,:) + cYv;
%                             dY(iNode,jNode,:) = dY(iNode,jNode,:) + dYv;
                    %%%%

                            

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