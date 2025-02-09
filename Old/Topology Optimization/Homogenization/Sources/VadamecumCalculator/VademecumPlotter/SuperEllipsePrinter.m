classdef SuperEllipsePrinter < handle
    
    properties (Access = private)
        mesh
        outPutName
        iteration
    end
    
    methods (Access = public)
        
        function obj = SuperEllipsePrinter(cParams)
            obj.init(cParams)            
        end
        
        function print(obj)
            dI.mesh            = obj.mesh;
            dI.outName         = obj.outPutName;
            dI.pdim            = '2D';
            dI.ptype           = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
            dB.printers = 'Density';
            postCase = 'Density';
            postProcess = Postprocess(postCase,dB);
            d.x = ones(size(obj.mesh.coord(:,1),1),1);
            postProcess.print(obj.iteration,d);
        end
        
        function captureImage(obj)
            f = obj.outPutName;
            it = obj.iteration;
            imagePath = '/home/alex/git-repos/MicroStructurePaper/';
            outPutNameImage = fullfile(imagePath,[f,num2str(it)]);
            inputFileName = fullfile('Output',f,[f,num2str(it),'.flavia.res']);
            sI.fileName = f;
            sI.outPutImageName = outPutNameImage;
            sI.inputFileName = inputFileName;
            imageCapturer = GiDImageCapturer(sI);
            imageCapturer.capture();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.outPutName = cParams.outPutName;
            obj.mesh       = cParams.mesh;   
            obj.iteration  = 0;
        end
        
    end
    
end