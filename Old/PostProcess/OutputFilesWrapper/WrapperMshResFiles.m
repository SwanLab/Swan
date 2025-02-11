classdef WrapperMshResFiles < handle
    
    properties (Access = public)
        dataRes
        mesh
    end
    
    properties (Access = private)
        dataMesh
    end
    
    properties (Access = private)
        fileName
        folderPath
    end
    
    methods (Access = public)
        
        function obj = WrapperMshResFiles(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj)
            obj.wrapGiDMshData();
            obj.createMesh();
            obj.wrapGiDResData();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName   = cParams.fileName;
            obj.folderPath = cParams.folderPath;
        end
        
        function wrapGiDMshData(obj)
            fM = [obj.fileName,'.flavia.msh'];
            s.filePath = fullfile(obj.folderPath,fM);
            wM = WrapperMshFile(s);
            wM.read();
            obj.dataMesh = wM.getDataBase();
        end
        
        function createMesh(obj)
            s.connec = obj.dataMesh.connec;
            s.coord  = obj.dataMesh.coord;
            obj.mesh = Mesh.create(s);
        end
        
        function wrapGiDResData(obj)
            fR = [obj.fileName,'.flavia.res'];
            s.filePath = fullfile(obj.folderPath,fR);                        
            s.nElem     = obj.mesh.nelem;
            s.nNodes    = obj.mesh.nnodes;
            s.dimension = obj.mesh.ndim;
            wR = WrapperResFile(s);
            wR.read();
            obj.dataRes = wR.getDataBase();
        end
        
    end

    
end