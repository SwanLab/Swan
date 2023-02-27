classdef TotalCorrectorComputer < handle
    
    
    properties (Access = private)
        ortoghonalCorrectors  
        ortogonalCoefficients
        totalCorrector
    end
    
    properties (Access = private)
       mesh 
       orientationVector
       interpolator
       dilatedOrientation
       singularities
       phiMapping
       isCoherent
    end
    
    methods (Access = public)

        function obj = TotalCorrectorComputer(cParams)
            obj.init(cParams); 
            obj.createTotalCorrector();
        end

        function computeCoefficients(obj)
            nSing = size(obj.singularities,1)-1;
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
                for iSing = 1:size(obj.singularities)-1
                    ocV = oC{iSing}.fValues;
                    psiT(iDim,:,:) = psiT(iDim,:,:) + coef(iDim,iSing)*ocV;
                end
            end
            obj.totalCorrector.fValues = psiT;
        end

        function tC = compute(obj)
           if ~isempty(obj.singularities)               
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
           obj.singularities      = cParams.singularities;
           obj.isCoherent         = cParams.isCoherent;
           obj.interpolator       = cParams.interpolator;
           obj.dilatedOrientation = cParams.dilatedOrientation;
           obj.phiMapping         = cParams.phiMapping;
        end

        function createTotalCorrector(obj)
            nElem = obj.mesh.nelem;
            ndif  = obj.mesh.ndim;
            nnode = obj.mesh.nnodeElem;
            s.fValues = zeros(ndif,nnode,nElem);
            s.mesh    = obj.mesh;
            obj.totalCorrector = P1DiscontinuousFunction(s);    
        end                   
        
        function computeOrtoghonalCorrector(obj)
            nSing = size(obj.singularities,1)-1;
            oC = cell(nSing,1);
            for iS = 1:nSing
                sCoord = obj.singularities(iS,:);
                cF     = obj.computeCorrectorFunction(sCoord);
                sF     = obj.createShifting(cF);
                oC{iS} = obj.computeOrthogonalCorrector(cF,sF);
            end
            obj.ortoghonalCorrectors = oC;
        end                

        function cV = computeCorrectorFunction(obj,sCoord)
            s.mesh             = obj.mesh;
            s.isCoherent       = obj.isCoherent;
            s.singularityCoord = sCoord;
            c = CorrectorComputer(s);
            cV = c.compute();
        end            
        
        function sF = createShifting(obj,cF)
            s.mesh         = obj.mesh;
            s.corrector    = cF;
            s.interpolator = obj.interpolator;
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