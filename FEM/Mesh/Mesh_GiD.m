classdef Mesh_GiD < Mesh

    properties (GetAccess = public,SetAccess = public)
        pdim
        dirichlet
        pointload
    end
    
    properties (GetAccess = public,SetAccess = public)
        ptype
        scale
        problemID
    end
    
    properties (Access = public)
        femReader
    end
    
    methods (Access = public)
        
        function obj = Mesh_GiD(fileName)
            obj.readFemInput(fileName);
            obj.assignParams();
            obj.createMesh();
        end
        
        function copy = clone(obj)
            copy = Mesh_GiD(obj.problemID);
        end
        
        function simpleCopy = getSimplifiedMesh(obj)
            simpleCopy = Mesh();
            simpleCopy.create(obj.coord,obj.connec);
        end
        
    end
    
    methods (Access = public)
        
        function readFemInput(obj,fileName)
            obj.femReader = FemInputReader_GiD();
            obj.femReader.read(fileName);
        end
        
        function assignParams(obj)
            obj.pdim = obj.femReader.pdim;
            obj.dirichlet = obj.femReader.dirichlet;
            obj.pointload = obj.femReader.pointload;
            obj.ptype = obj.femReader.ptype;
            obj.scale = obj.femReader.scale;            
            obj.problemID = obj.femReader.problemID;
            obj.geometryType = obj.femReader.geometryType;
            obj.coord = obj.femReader.coord;            
            obj.connec = obj.femReader.connec;
        end
        
        function createMesh(obj)
            obj.create(obj.coord,obj.connec);
        end
        
    end 
    
end

