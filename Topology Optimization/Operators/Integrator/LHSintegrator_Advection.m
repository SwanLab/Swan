classdef LHSintegrator_Advection < handle

    properties (Access = private)
        mesh
        test, trial
        b
        rho
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegrator_Advection(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = compute(obj)
            [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = obj.computeElementalLHS();
            CX = obj.assembleMatrix(CX);
            CY = obj.assembleMatrix(CY);
            DX = obj.assembleMatrix(DX);
            DY = obj.assembleMatrix(DY);
            EX = obj.assembleMatrix(EX);
            EY = obj.assembleMatrix(EY);
            Kxx = obj.assembleMatrix(Kxx);
            Kxy = obj.assembleMatrix(Kxy);
            Kyy = obj.assembleMatrix(Kyy);
            Mxx = obj.assembleMatrix(Mxx);
            Mxy = obj.assembleMatrix(Mxy);
            Myy = obj.assembleMatrix(Myy);
            Mrho = obj.assembleMatrix(Mrho);
        end

    end

    methods (Access = protected)

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            dNdxTs = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodETs = size(dNdxTs,2);
            nDofETs = nNodETs*obj.test.ndimf;
            nNodETr = size(dNdxTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;

            shapesTest  = obj.test.computeShapeFunctions(xV);
            shapesTrial = obj.trial.computeShapeFunctions(xV);

            % N  = obj.interpolation.shape;


            bD   = obj.b.project('P1D');
            rhoD = obj.rho.project('P1D');

            bElem   = bD.fValues;
            rhoElem = rhoD.fValues;



            cX = zeros(nDofETs,nDofETr,nElem);
            cY = zeros(nDofETs,nDofETr,nElem);
            dX = zeros(nDofETs,nDofETr,nElem);
            dY = zeros(nDofETs,nDofETr,nElem);
            eX = zeros(nDofETs,nDofETr,nElem);
            eY = zeros(nDofETs,nDofETr,nElem);
            Kxx = zeros(nDofETs,nDofETr,nElem);
            Kyy = zeros(nDofETs,nDofETr,nElem);
            Kxy = zeros(nDofETs,nDofETr,nElem);
            Mxx = zeros(nDofETs,nDofETr,nElem);
            Myy = zeros(nDofETs,nDofETr,nElem);
            Mxy = zeros(nDofETs,nDofETr,nElem);
            Mrho = zeros(nDofETs,nDofETr,nElem);

            for igaus = 1:nGaus
                dNg     = dNdxTs(:,:,:,igaus);
                Ng      = shapesTest(:,igaus);
                dV(:,1) = dVolu(igaus, :);
                for iNode = 1:nNodETs
                    dNi = squeezeParticular(dNg(:,iNode,:),2);
                    Ni  = Ng(iNode);
                    for jNode = 1:nNodETs
                        dNj = squeezeParticular(dNg(:,jNode,:),2);
                        Nj  = Ng(jNode);
                        for kNode = 1:nNodETs
                            bXk = squeeze(bElem(1,kNode,:));
                            bYk = squeeze(bElem(2,kNode,:));
                            rhok = squeeze(rhoElem(1,kNode,:));

                            Nk  = Ng(kNode);
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

                            for lNode = 1:nNodETs


                                Nl = Ng(lNode);
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
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.b            = cParams.b;
            obj.rho          = cParams.rho;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.trial.order;
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assembleFunctions(lhs, obj.test, obj.trial);
        end

    end

end