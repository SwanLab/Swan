classdef computingOptimal < handle
    
    methods (Access = public)
        
        function obj = computingOptimal()
            obj.init()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function obj = computingOpt()
            fileName = 'OptimalSuperEllipse';
            outputFolder = fullfile(pwd,'Output',fileName);
            gmsFile = [fullfile(outputFolder,fileName),'.msh'];
            createMesh(fileName,outputFolder);
            homog = obj.createNumericalHomogenizer(fileName,gmsFile);
            maxSnorm = obj.computeMaxStressNorm(homog);
        end
        
        
        
        function maxSnorm = computeMaxStressNorm(obj,homog)
            Ch = homog.cellVariables.Ch;
            microProblem = homog.getMicroProblem();
            stress = [cos(pi/4) sin(pi/4) 0]';
            strain = Ch\stress;
            microProblem.element.setVstrain(strain');
            microProblem.computeVariables;
            stresses = microProblem.variables.stress;
            sx = squeeze(stresses(:,1,:));
            sy = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            maxSnorm = max(sNorm);
        end
        
        function createMesh(obj,fileName,outputFolder)
            d = SettingsFreeFemMeshGenerator();
            d.freeFemFileName = 'SmoothRectangle';
            d.hMax  = 0.02;%0.0025;
            d.mxV             = 0.5;
            d.myV             = 0.5;
            d.fileName        = fileName;
            d.printingDir     = outputFolder;
            d.qNorm           = 32;
            
            fG = FreeFemMeshGenerator(d);
            fG.generate();
            
        end
        
        function homog = createNumericalHomogenizer(obj,fileName,gmsFile)
            d.gmsFile = gmsFile;
            d.outFile = fileName;
            d.print   = false;
            d.iter = 0;
            nH = NumericalHomogenizerCreatorFromGmsFile(d);
            homog = nH.getHomogenizer();
        end
        
    end
    
end



