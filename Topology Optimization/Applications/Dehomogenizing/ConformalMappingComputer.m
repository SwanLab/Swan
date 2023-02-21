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
            for iDim = 1:obj.mesh.ndim
                obj.phi{iDim}.plotContour();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.orientation = cParams.theta;
            obj.mesh        = cParams.mesh;
            obj.dilation    = cParams.dilation;
        end

        function createInterpolator(obj)
            s.mesh        = obj.mesh;
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
            %nnod = m.nnodes;
            nnod = obj.mesh.nelem*obj.mesh.nnodeElem;   
           % phiI = zeros(nnod,nDim);

            if ~isempty(obj.singularityCoord)
                for iS = 1:size(obj.singularityCoord)-1
                    sCoord = obj.singularityCoord(iS,:);
                    b1 = obj.orientationVector{1}.fValues;
                    cr = obj.computeCorrector(b1,sCoord);
                    oC{iS} = obj.computeOrthogonalCorrector(cr);
                end
            end
    
            
            for iDim = 1:nDim
                bI    = obj.orientationVector{iDim}.fValues;
                phiD  = obj.computeMapping(bI);
%                phiI(:,iDim) = phiD.getFvaluesAsVector();
                phiI(iDim,:,:) = phiD.fValues;
            end

             if ~isempty(obj.singularityCoord)
                    psiT = zeros(size(phiI));
                 for iDim = 1:nDim            
                    bI = obj.orientationVector{iDim}.fValues;
                    coef = obj.computeCoeffs(oC,bI);
                    
                    for iSing = 1:size(obj.singularityCoord)-1
                        ocV = oC{iSing}.fValues;                        
                        psiT(iDim,:,:) = psiT(iDim,:,:) + coef(iSing)*ocV;
                    end
                 end 
                 phiI = phiI + psiT;                 
            end 

            s.fValues = phiI;
            s.mesh    = obj.mesh;
            psiTs = P1DiscontinuousFunction(s);
            obj.phi = psiTs;
        end

        function phi = computeMapping(obj,fValues)
            s.fValues = fValues';
            s.mesh    = obj.mesh;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            problem   = MinimumDiscGradFieldWithVectorInL2(s);
            phi = problem.solve();
        end

        function cV = computeCorrector(obj,b,sCoord)
            s.mesh               = obj.mesh;
            s.orientation        = b;
            s.singularityCoord = sCoord;
            c = CorrectorComputer(s);
            cV = c.compute();
            % c.plot()
        end


        function computeOrientationVector(obj)
            er = exp(obj.dilation);
            erCos = er.*cos(obj.orientation);
            erSin = er.*sin(obj.orientation);
            b1(:,1) = erCos;
            b1(:,2) = erSin;
            b2(:,1) = -erSin;
            b2(:,2) = erCos;
            b(:,:,1) = b1;
            b(:,:,2) = b2;
            for iDim = 1:obj.mesh.ndim
                s.fValues = b(:,:,iDim);
                s.mesh   = obj.mesh;
                bf = P1Function(s);
                obj.orientationVector{iDim} = bf;
            end
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

            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');

            nSing = numel(psi);
            for iSing = 1:nSing
                dPsiV = psi{iSing}.computeGradient(q);
                dPsi(:,iSing,:,:) = dPsiV.fValues;
            end

            LHS = zeros(nSing,nSing);
            RHS = zeros(nSing,1);


            s.fValues = b;
            s.mesh   = obj.mesh;
            bf = P1Function(s);


            xGauss = q.posgp;
            bfG    = bf.evaluate(xGauss);



            dV = obj.mesh.computeDvolume(q);
            dVt(1,:,:) = dV;
            dVT = repmat(dVt,bf.ndimf,1,1);

            for iSing = 1:nSing
                dPsiI = permute(squeeze(dPsi(:,iSing,:,:)),[1 3 2]);
                for jSing = 1:nSing
                    dPsiJ = permute(squeeze(dPsi(:,jSing,:,:)),[1 3 2]);
                    lhs   = dPsiI.*dPsiJ.*dVT;
                    LHS(iSing,jSing) = sum(lhs(:));
                end
                rhs = dPsiI.*bfG.*dVT;
                RHS(iSing) = sum(rhs(:));
            end
           c = LHS\RHS;
        end

    end

end
