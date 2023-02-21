classdef TotalCorrectorComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        ortoghonalCorrector  
        ortogonalCoefficients
        totalCorrector
    end
    
    properties (Access = private)
       mesh 
       orientationVector
       singularityCoord
       interpolator
       phiMapping
    end
    
    methods (Access = public)

        function obj = TotalCorrectorComputer(cParams)
            obj.init(cParams); 
            obj.createTotalCorrector();
        end

        function computeCoefficients(obj)
            nSing = size(obj.singularityCoord,1)-1;
            nDim  = obj.mesh.ndim;
            c = zeros(nDim,nSing);
            for iDim = 1:nDim
                bI = obj.orientationVector{iDim}.fValues;
                oC = obj.ortoghonalCorrector;
                c(iDim,:) = obj.computeCoeffs(oC,bI);
            end
            obj.ortogonalCoefficients = c;
        end

        function psiT = computeTotalComputer(obj)
            oC   = obj.ortoghonalCorrector;
            coef = obj.ortogonalCoefficients;
            psiT = obj.totalCorrector.fValues;
            for iDim = 1:obj.mesh.ndim
                for iSing = 1:size(obj.singularityCoord)-1
                    ocV = oC{iSing}.fValues;
                    psiT(iDim,:,:) = psiT(iDim,:,:) + coef(iDim,iSing)*ocV;
                end
            end
            obj.totalCorrector.fValues = psiT;
        end

        function tC = compute(obj)
           obj.computeSingularities();           
           if ~isempty(obj.singularityCoord)               
               obj.computeOrtoghonalCorrector();
               obj.computeCoefficients();
               obj.computeTotalComputer();
           end    
           tC = obj.totalCorrector;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh = cParams.mesh;
           obj.orientationVector = cParams.orientationVector;
           obj.interpolator = cParams.interpolator;
           obj.phiMapping   = cParams.phiMapping;
        end

        function createTotalCorrector(obj)
            s.fValues = zeros(size(obj.phiMapping.fValues));
            s.mesh    = obj.mesh;
            obj.totalCorrector = P1DiscontinuousFunction(s);    
        end       
        

        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientationVector{1};
            sC = SingularitiesComputer(s);
            sCoord = sC.compute();
            %sC.plot();
            obj.singularityCoord = sCoord;
        end        
        
        function computeOrtoghonalCorrector(obj)
            nSing = size(obj.singularityCoord,1)-1;
            oC = cell(nSing,1);
            for iS = 1:nSing
                sCoord = obj.singularityCoord(iS,:);
                b1 = obj.orientationVector{1}.fValues;
                cr = obj.computeCorrector(b1,sCoord);
                oC{iS} = obj.computeOrthogonalCorrector(cr);
            end
            obj.ortoghonalCorrector = oC;
        end        
        
        function c = computeOrthogonalCorrector(obj,c)
            s.mesh               = obj.mesh;
            s.correctorValue     = c;
            s.interpolator       = obj.interpolator;
            o = OrthogonalCorrectorComputer(s);
            c = o.compute();
            % o.plot();
        end    
        

        function cV = computeCorrector(obj,b,sCoord)
            s.mesh               = obj.mesh;
            s.orientation        = b;
            s.singularityCoord = sCoord;
            c = CorrectorComputer(s);
            cV = c.compute();
            % c.plot()
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