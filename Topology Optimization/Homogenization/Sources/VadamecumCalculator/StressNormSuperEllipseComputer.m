classdef StressNormSuperEllipseComputer < handle
        
    properties (Access = private)
        homog        
        fName        
        outputFolder                
    end
    
    properties (Access = private)
        mx
        my
        q
        phi
        pNorm
        print
        hMesh
        fileName
        hasToCaptureImage
    end
    
    methods (Access = public)
        
        function obj = StressNormSuperEllipseComputer(cParams)
            obj.init(cParams)
        end        
        
        function sPnorm = compute(obj)
            obj.createMesh();
            obj.createNumericalHomogenizer();
            sPnorm = obj.computePstressNorm();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mx       = cParams.mx;
            obj.my       = cParams.my;
            obj.q        = cParams.q;
            obj.phi      = cParams.phi;
            obj.pNorm    = cParams.pNorm;
            obj.print    = cParams.print;
            obj.hMesh    = cParams.hMesh; 
            obj.fileName = cParams.fileName;
            obj.hasToCaptureImage = cParams.hasToCaptureImage;
            obj.computeFileNameAndOutputFolder();
        end
        
        function computeFileNameAndOutputFolder(obj)
            obj.fName = [obj.fileName];
            obj.outputFolder = fullfile(pwd,'Output',obj.fName);
        end
        
        function sPnorm = computePstressNorm(obj)
            Ch = obj.homog.cellVariables.Ch;
            microProblem = obj.homog.getMicroProblem();
            stress = [cos(obj.phi) sin(obj.phi) 0]';
            strain = Ch\stress;
            microProblem.element.setVstrain(strain');
            microProblem.computeVariables;
            stresses = microProblem.variables.stress;
            m = microProblem.mesh;
            quad = Quadrature.set(m.geometryType);
            quad.computeQuadrature('CONSTANT');
            dV = m.computeDvolume(quad);
            sx  = squeeze(stresses(:,1,:));
            sy  = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm2 = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            if isequal(obj.pNorm,'max')
                sPnorm = max(sNorm2);
            else
                p = obj.pNorm;
                int = sNorm2.^p;
                sPnorm = sum(int(:).*dV(:))^(1/p);
            end
        end
        
        function createNumericalHomogenizer(obj)
            d.gmsFile = [fullfile(obj.outputFolder,obj.fName),'.msh'];
            d.outFile = obj.fName;
            d.print   = obj.print;
            d.iter = 0;
            d.hasToCaptureImage = obj.hasToCaptureImage;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            obj.homog = nH.getHomogenizer();
        end
                
        function createMesh(obj)
            d = SettingsFreeFemMeshGenerator();
            d.freeFemFileName = 'SmoothRectangle';
            d.hMax  = obj.hMesh;%0.002;%0.0025;
            d.mxV             = obj.mx;
            d.myV             = obj.my;
            d.fileName        = obj.fileName;
            d.printingDir     = obj.outputFolder;
            d.qNorm           = obj.q;
            fG = FreeFemMeshGenerator(d);
            fG.generate();
        end        
        
    end
    
        
end