classdef DehomogenisationMeshPrinter < FilePrinter
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = DehomogenisationMeshPrinter(cParams)
            obj.init(cParams)            
        end
        
        function print(obj)
            obj.openFile();
            obj.printLines();
            obj.closeFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.mesh     = cParams.mesh;
        end
        
        function printLines(obj)
            obj.printCoord();
            obj.printConnec();
        end
        
        function printCoord(obj)
           fprintf(obj.fileID,'nodes coordinates \n'); 
           formatV = repmat('%4.12f ',1,size(obj.mesh.coord,2));
           format = [formatV,'\n'];
           fprintf(obj.fileID,format,obj.mesh.coord');            
           fprintf(obj.fileID,'\n'); 
        end
        
        function printConnec(obj)
           fprintf(obj.fileID,'connectivities \n'); 
           formatV = repmat('%12.0f ',1,size(obj.mesh.connec,2));
           format = [formatV,'\n'];
           fprintf(obj.fileID,format,obj.mesh.connec');
        end       
        
    end
    
end