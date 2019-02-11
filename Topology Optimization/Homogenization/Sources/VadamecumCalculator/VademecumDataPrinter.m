classdef VademecumDataPrinter < FilePrinter
    
    properties (Access = private)
        Ctensor
        Ptensor
        tensor
        tensorName
        tensComp
        outPutPath
    end
    
    
    methods (Access = public)
        
        
        function obj = VademecumDataPrinter(d)
            obj.init(d);
        end
        
        function print(obj)
            obj.printCtensor();
            obj.printPtensor();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.outPutPath = d.outPutPath;
            obj.Ctensor    = d.Ctensor;
            obj.Ptensor    = d.Ptensor;
        end
        
        function printCtensor(obj)
            obj.tensor     = obj.Ctensor;
            obj.tensorName = 'C_' ;
            obj.printTensor();
        end
        
        function printPtensor(obj)
            obj.tensor     = obj.Ptensor;
            obj.tensorName = 'P_' ;
            obj.printTensor();
        end
        
        function printTensor(obj)
            [nx,ny] = size(obj.tensor(:,:,1,1));
            for i = 1:nx
                for j = 1:ny
                    obj.obtainTensorComponent(i,j);
                    obj.obtainFileName(i,j);
                    obj.openFile();
                    obj.printTensorComponent();
                    obj.closeFile();
                end
            end
        end
        
        function obtainTensorComponent(obj,i,j)
            T = obj.tensor(i,j,:,:);
            obj.tensComp = squeeze(T);
        end
        
        function obtainFileName(obj,i,j)
            path = obj.outPutPath;
            tName = obj.tensorName;
            istr = num2str(i);
            jstr = num2str(j);
            obj.fileName = [path,tName,istr,jstr,'.txt'];
        end
        
        function printTensorComponent(obj)
            T  = obj.tensComp;
            fN = obj.fileName;
            save(fN,'T','-ascii','-double','-tabs');
        end
        
        
    end
    
end