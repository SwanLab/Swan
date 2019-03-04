classdef VademecumDataPrinter < FilePrinter
    
    properties (Access = private)
        Ctensor
        Ptensor
        volume
        tensor
        fieldName
        fieldValue
        outPutPath
    end
    
    
    methods (Access = public)
        
        
        function obj = VademecumDataPrinter(d)
            obj.init(d);
        end
        
        function print(obj)
            obj.printVolume();
            obj.printCtensor();
            obj.printPtensor();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.outPutPath = d.outPutPath;
            obj.Ctensor    = d.Ctensor;
            obj.Ptensor    = d.Ptensor;
            obj.volume     = d.volume;
        end
        
        function printVolume(obj)
            obj.fieldValue = obj.volume;
            tName = 'Volume' ;
            path  = obj.outPutPath;
            obj.fileName = [path,tName,'.txt'];
            obj.openFile();
            obj.printVariable();
            obj.closeFile();            
        end
        
        function printCtensor(obj)
            obj.tensor     = obj.Ctensor;
            obj.fieldName = 'C_' ;
            obj.printTensor();
        end
        
        function printPtensor(obj)
            obj.tensor     = obj.Ptensor;
            obj.fieldName = 'P_' ;
            obj.printTensor();
        end
        
        function printTensor(obj)
            [nx,ny] = size(obj.tensor(:,:,1,1));
            for i = 1:nx
                for j = 1:ny
                    obj.obtainTensorComponent(i,j);
                    obj.obtainTensorFileName(i,j);
                    obj.openFile();
                    obj.printVariable();
                    obj.closeFile();
                end
            end
        end
        
        function obtainTensorComponent(obj,i,j)
            t = obj.tensor(i,j,:,:);
            obj.fieldValue = squeeze(t);
        end
        
        function obtainTensorFileName(obj,i,j)
            path  = obj.outPutPath;
            tName = obj.fieldName;
            istr = num2str(i);
            jstr = num2str(j);
            obj.fileName = [path,tName,istr,jstr,'.txt'];
        end
        
        function printVariable(obj)
            fV  = obj.fieldValue;
            fN  = obj.fileName;
            save(fN,'fV','-ascii','-double','-tabs');
        end
        
        
    end
    
end