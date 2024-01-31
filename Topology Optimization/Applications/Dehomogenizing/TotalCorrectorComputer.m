classdef TotalCorrectorComputer < handle
    
    
    properties (Access = private)
        nCorr        
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
            obj.computeNumberOfCorrectors();
        end

        function tC = compute(obj)
           if obj.nCorr > 0               
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

        function computeNumberOfCorrectors(obj)
            obj.nCorr =  obj.singularities.nSing - 1;
        end

        function computeOrtoghonalCorrector(obj)
            oC      = cell(obj.nCorr,1);
            areSing = find(squeeze(obj.singularities.isElemSingular.fValues));           
            for iS = 1:obj.nCorr
              %  sCoord = obj.singularities.coord(iS,:);
           %   iSing  = areSing(end+1-iS);
                iSing  = areSing(iS);
                cF     = obj.computeCorrectorFunction(iSing);
                sF     = obj.createShifting(cF);
                oC{iS} = obj.computeOrthogonalCorrector(cF,sF);
            end
            obj.ortoghonalCorrectors = oC;
        end             
        
       function computeCoefficients(obj)
            nDim  = obj.mesh.ndim;
            c = zeros(nDim,obj.nCorr);
            for iDim = 1:nDim
                bI = obj.dilatedOrientation{iDim};
                oC = obj.ortoghonalCorrectors;
                c(iDim,:) = obj.computeCoeffs(oC,bI);
            end
            obj.ortogonalCoefficients = c;
        end

        function psiT = computeTotalComputer(obj)
            oC    = obj.ortoghonalCorrectors;
            coef  = obj.ortogonalCoefficients;
            psiT  = obj.totalCorrector.fValues;
            for iDim = 1:obj.mesh.ndim
                for iSing = 1:obj.nCorr
                    ocV = oC{iSing}.fValues;
                    psiT(iDim,:,:) = psiT(iDim,:,:) + coef(iDim,iSing)*ocV;
                end
            end
            obj.totalCorrector.fValues = psiT;
        end        

        function cV = computeCorrectorFunction(obj,singElem)
            s.mesh               = obj.mesh;
            s.isCoherent         = obj.isCoherent;
            s.singularElement    = singElem;
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
