classdef TotalCorrectorComputer < handle
    
    
    properties (Access = private)
        ortoghonalCorrectors  
        ortogonalCoefficients
        totalCorrector
    end
    
    properties (Access = private)
       mesh 
       dilatedOrientation
       orientationVector
       singularityCoord
       interpolator
       phiMapping
       isCoherent
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
                bI = obj.dilatedOrientation{iDim};
                oC = obj.ortoghonalCorrectors;
                c(iDim,:) = obj.computeCoeffs(oC,bI);
            end
            obj.ortogonalCoefficients = c;
        end

        function psiT = computeTotalComputer(obj)
            oC   = obj.ortoghonalCorrectors;
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
           obj.mesh               = cParams.mesh;
           obj.orientationVector  = cParams.orientationVector;
           obj.dilatedOrientation = cParams.dilatedOrientation;
           obj.phiMapping         = cParams.phiMapping;
        end

        function createTotalCorrector(obj)
            s.fValues = zeros(size(obj.phiMapping.fValues));
            s.mesh    = obj.mesh;
            obj.totalCorrector = P1DiscontinuousFunction(s);    
        end               

        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientationVector.value{1};
            sC = SingularitiesComputer(s);
            sCoord = sC.compute();
            obj.singularityCoord = sCoord;
        end        
        
        function computeOrtoghonalCorrector(obj)
            nSing = size(obj.singularityCoord,1)-1;
            oC = cell(nSing,1);
            for iS = 1:nSing
                sCoord = obj.singularityCoord(iS,:);
                cF     = obj.computeCorrectorFunction(sCoord);
                sF     = obj.createShifting(cF);
                oC{iS} = obj.computeOrthogonalCorrector(cF,sF);
            end
            obj.ortoghonalCorrectors = oC;
        end                

        function cV = computeCorrectorFunction(obj,sCoord)
            s.mesh              = obj.mesh;
            s.orientationVector = obj.orientationVector;
            s.singularityCoord = sCoord;
            c = CorrectorComputer(s);
            cV = c.compute();
        end            
        
        function sF = createShifting(obj,cF)
            s.mesh         = obj.mesh;
            s.corrector    = cF;
            s.interpolator = obj.orientationVector.interpolator;
            m = ShiftingFunctionComputer(s);
            sF = m.compute();
        end   

        function oC = computeOrthogonalCorrector(obj,cF,sF)
            phi = cF.fValues;
            fD  = sF.fValues;
            phi = phi - fD;
            s.mesh    = obj.mesh;
            s.fValues = phi;
            oC = P1DiscontinuousFunction(s);  
        end          
        
        function coef = computeCoeffs(obj,psi,b)
            s.orthogonalCorrector = psi;
            s.mesh                = obj.mesh;
            c = CorrectorCoefficientsComputer(s);
            coef = c.compute(b);
        end
        
    end
    
end