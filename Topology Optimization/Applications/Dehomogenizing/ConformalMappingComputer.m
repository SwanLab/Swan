classdef ConformalMappingComputer < handle


    properties (Access = private)
       phi
       interpolator
       singularityCoord
    end

    properties (Access = private)
       orientation
       orientationVector
       mesh
       dilation
    end

    methods (Access = public)

        function obj = ConformalMappingComputer(cParams)
            obj.init(cParams);
            obj.computeOrientationVector();
            obj.createInterpolator();
        end

        function phiV = compute(obj)
            obj.computeSingularities();
            obj.computeMappingWithSingularities();
            phiV = obj.phi;
        end

        function plot(obj)
           phi1 = obj.phi(:,1);
           phi2 = obj.phi(:,2);
           obj.plotContour((phi1));
           obj.plotContour((phi2));
        end

    end

    methods (Access = private)

        function plotContour(obj,z)
            figure()
            m = obj.mesh.createDiscontinousMesh;
            x = m.coord(:,1);
            y = m.coord(:,2);
            [~,h] = tricontour(m.connec,x,y,z,30);
            set(h,'LineWidth',5);
            colorbar
        end

        function init(obj,cParams)
            obj.orientation = cParams.theta;
            obj.mesh        = cParams.mesh;
            obj.dilation    = cParams.dilation;
        end

        function createInterpolator(obj)
            s.meshCont    = obj.mesh;
            s.meshDisc    = obj.mesh.createDiscontinuousMesh();
            s.orientation = [cos(obj.orientation),sin(obj.orientation)];
            s = SymmetricContMapCondition(s);
            sC = s.computeCondition();
            obj.interpolator = sC;
        end

        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = [cos(obj.orientation),sin(obj.orientation)];%obj.computeOrientationVector(1)';
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
           % sF.plot();
            coordB = obj.mesh.computeBaricenter();
            coordB = transpose(coordB);
            sCoord =  coordB(isS,:);
            obj.singularityCoord = sCoord;
        end

        function computeMappingWithSingularities(obj)
           nDim = 2;
           m = obj.mesh.createDiscontinuousMesh();
           nnod = m.nnodes;
           phiI = zeros(nnod,nDim);


           for iDim = 1:nDim
                b  = obj.orientationVector;
                bI = squeezeParticular(b(:,iDim,:),2);
                phiD  = obj.computeMapping(bI);     
                phiI(:,iDim) = phiD.getFvaluesAsVector();
                %phiI(:,iDim) = reshape(phiD',[],1);
           end

           if ~isempty(obj.singularityCoord)
               for iS = 1:size(obj.singularityCoord)-1
                sCoord = obj.singularityCoord(iS,:);
                b = obj.orientationVector;
                b1 = squeezeParticular(b(:,1,:),2);
                cr = obj.computeCorrector(b1,sCoord);
                oC{iS} = obj.computeOrthogonalCorrector(cr)
               end
           end         

           psiTs = zeros(nnod,nDim);
           if ~isempty(obj.singularityCoord)
               for iDim = 1:nDim
                   b  = obj.orientationVector;
                   bI = squeezeParticular(b(:,iDim,:),2);
                   coef = obj.computeCoeffs(oC,bI);
                   % coef = floor(coef);
                   psiT = zeros(size(phiI,1),1);
                   for iSing = 1:size(obj.singularityCoord)-1
                       ocV = oC{iSing}.getFvaluesAsVector();
                       psiT = psiT + coef(iSing)*ocV;
                   end
                   psiTs(:,iDim) = psiT;
               end
           end
           phiIc2  = phiI + psiTs;
           obj.phi = phiIc2;
        end

        function phi = computeMapping(obj,fValues)
            s.fValues = fValues;
            s.mesh    = obj.mesh;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            problem   = MinimumDiscGradFieldWithVectorInL2(s);
            phi = problem.solve();
        end

        function cV = computeCorrector(obj,b,sCoord)            
            s.mesh               = obj.mesh;
            s.orientation        = b';
            s.singularityCoord = sCoord;
            c = CorrectorComputer(s);
            cV = c.compute();
           % c.plot()
        end


        function computeOrientationVector(obj)
           er = exp(obj.dilation);
           erCos = er.*cos(obj.orientation);
           erSin = er.*sin(obj.orientation);
           Q(1,1,:) = erCos;
           Q(1,2,:) = -erSin;
           Q(2,1,:) = erSin;
           Q(2,2,:) = erCos;
           obj.orientationVector = Q;
        end

        function c = computeOrthogonalCorrector(obj,c)
            s.mesh               = obj.mesh;
            s.correctorValue     = c;
            s.interpolator       = obj.interpolator;
            o = OrthogonalCorrectorComputer(s);
            c = o.compute();
           % o.plot();
        end

        function c = computeCoeffs(obj,psi,b)
            bGauss = obj.computeFGauss(b);
            I = obj.interpolator;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            m = obj.mesh.createDiscontinuousMesh();
            dV = obj.mesh.computeDvolume(q);
            nDim = m.ndim;
            nSing = numel(psi);
            LHS = zeros(nSing,nSing);
            RHS = zeros(nSing,1);

            nnode = obj.mesh.nnodeElem;
            nElem = obj.mesh.nelem;
            dPsi = zeros(nDim,nnode,nElem,nSing);
            for iSing = 1:nSing
                psiV = squeeze(psi{iSing}.fValues)';
                dPsi(:,:,:,iSing) = obj.createDiscontinousGradient(psiV);
            end


            for igaus = 1:q.ngaus
                dVG  = dV(igaus,:)';
                for idim = 1:nDim
                    bG = squeeze(bGauss(idim,igaus,:));
                    dPsiG = squeeze(dPsi(idim,igaus,:,:));
                    for iSing = 1:nSing
                        dPsiI = dPsiG(:,iSing);
                        for jSing = 1:nSing
                            dPsiJ = dPsiG(:,:,:,jSing);
                            lhs = sum(dPsiI.*dVG.*dPsiJ);
                            LHS(iSing,jSing) = LHS(iSing,jSing) + lhs;
                        end
                        rhs = sum(dPsiI.*dVG.*bG);
                        RHS(iSing) = RHS(iSing) + rhs;
                    end

                end

            end

            LHS2 = zeros(nSing,nSing);
            RHS2 = zeros(nSing,1);


            sF.mesh               = m;
            sF.ndimf              = 1;
            sF.interpolationOrder = m.interpolation.order;
            s.field = Field(sF);
            s.mesh         = m;
            s.globalConnec = m.connec;
            s.type         = 'StiffnessMatrix';

            psiV = squeeze(psi{1}.fValues)';
            lhs = LHSintegrator.create(s);
            Kdg = lhs.compute();
            phidg = reshape(psiV',[],1);
            LHS2 = phidg'*Kdg*phidg;


            q = Quadrature.set(m.type);
            q.computeQuadrature('QUADRATIC');
            int = Interpolation.create(m,m.interpolation.order);
            int.computeShapeDeriv(q.posgp);
            s.mesh = m;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dN = g.dNdx;

            %dN = int.deriv
            inte = zeros(m.nelem,nnode);
            for idim = 1:nDim
                for igaus = 1:q.ngaus
                    dVG  = dV(igaus,:)';
                    bG = squeeze(bGauss(idim,igaus,:));
                    for iNode = 1:nnode
                        dNi = squeeze(dN(idim,iNode,:,igaus));
                        intG = dNi.*bG.*dVG;
                        inte(:,iNode) = inte(:,iNode) + intG;
                    end
                end
            end
            RHS2 = sum(sum(psiV.*inte));


            c = LHS2\RHS2;
        end

        function bG = computeFGauss(obj,b)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            xGauss = q.posgp;
            bG = zeros(obj.mesh.ndim,q.ngaus,obj.mesh.nelem);
            for idim = 1:obj.mesh.ndim
                s.fValues = b(idim,:)';
                s.connec = obj.mesh.connec;
                s.type   = obj.mesh.type;
                f = P1Function(s);
                bG(idim,:,:) = f.evaluate(xGauss);
            end
            %%%%% HEREEEE!!!!! Integrate with more gauss points b
        end        

        function fGauss = createDiscontinousGradient(obj,field)%%%Ehhhhh
            cV = field;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            int = Interpolation.create(obj.mesh,obj.mesh.interpolation.order);
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dN = g.dNdx;

            %dN = int.deriv;
            nDim = size(dN,1);
            nnode = size(dN,2);
            fGauss = zeros(nDim,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for idim = 1:nDim
                    for iNode = 1:nnode
                    dNi = squeeze(dN(idim,iNode,:,igaus));
            %        dNi = squeeze(dN(idim,iNode,igaus));

                    grad = dNi.*cV(:,iNode);
                    fG(:,1) = squeeze(fGauss(idim,igaus,:));
                    fGauss(idim,igaus,:) = fG + grad;
                    end
                end
            end
        end


    end

end
