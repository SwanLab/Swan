classdef NumericalRectangleTypeHomogenizer < NumericalHomogenizer
    
    properties (Access = protected, Abstract)
        outPutName
    end

    properties (Access = private)
        m1
        m2
    end
    
    methods (Access = public)
       
        function ls = getLevelSet(obj)
            ls = obj.levelSet;
        end
        
       function d = getDensity(obj)
            d = obj.density;                       
       end
       
       function f = getFileName(obj)
           f = obj.fileName;
       end
       
    end
    
    methods (Access = protected)
        
        function compute(obj,fileName,print,m1,m2,iter)
            obj.m1 = m1;
            obj.m2 = m2;
            obj.init(fileName,print,iter);
            obj.generateMicroProblem();
            obj.computeHomogenizedVariables();
            obj.createDensityPrinter()
            obj.print()
            obj.captureImage(print,iter)
        end
        
        function createDensity(obj)
            obj.createLevelSet(obj.m1,obj.m2)
            lsNodes = obj.levelSet;
            phyPr   = obj.microProblem;
            filter = FilterP0(lsNodes,phyPr);
            obj.density(:,1) = filter.getDens0();               
        end
        
    end
    
    methods (Access = private)
        
        function captureImage(obj,print,iter)
            if print
                f = obj.resFile;
                outPutNameWithIter = [obj.outPutName,num2str(iter)];
                inputFileName = fullfile('Output',f,[f,num2str(iter),'.flavia.res']);
                GiDImageCapturer(f,outPutNameWithIter,inputFileName);
            end
        end
        
    end
    
    methods (Access = protected, Abstract)
        createLevelSet(obj,m1,m2)
    end
    
end

