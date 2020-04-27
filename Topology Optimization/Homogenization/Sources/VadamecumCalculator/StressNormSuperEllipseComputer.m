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
        
        function printStress(obj)
            microProblem = obj.homog.getMicroProblem(); 
            dI.mesh    =  microProblem.mesh;
            dI.outName = [obj.fileName,'Print'];
            dI.pdim    = '2D';
            dI.ptype   = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
            postCase = 'ElasticityMicro';
            postProcess = Postprocess(postCase,dB);
            d.fields = microProblem.variables;
            d.quad = microProblem.element.quadrature;
            iter = 0;
            postProcess.print(iter,d);
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
            stress = [cos(obj.phi) sin(obj.phi) 0];
            microProblem = obj.homog.getMicroProblem();
            v = microProblem.computeVarFromStress(stress);
            stresses = v.stress;
            
            
            m = microProblem.mesh;
            quad = Quadrature.set(m.geometryType);
            quad.computeQuadrature('CONSTANT');
            dV = m.computeDvolume(quad);
            sx  = squeeze(stresses(:,1,:));
            sy  = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            if isequal(obj.pNorm,'max')
                sPnorm = max(sNorm);
            else
                p = obj.pNorm;
                int = sNorm.^p;
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