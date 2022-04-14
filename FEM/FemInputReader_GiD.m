classdef FemInputReader_GiD < handle
    
    properties (GetAccess = public,SetAccess = public)
        pdim
        dirichlet
        pointload
        ptype
        scale
        problemID
        fileName
    end
    
    properties (Access = public)
        geometryType
        coord
        connec
    end
    
    properties (Access = private)
        masterSlave
        boundaryNodes
        boundaryElements
    end
    
    methods (Access = public)
        
        function s = read(obj,fileName)
            obj.fileName = fileName;
            if ~isempty(fileName)
                obj.readFile(fileName);
            end
            s = obj.getData();
        end
        
        function s = getData(obj)
            s.mesh = obj.createMesh();
            s.pdim = obj.pdim;
            s.geometryType = obj.geometryType;
            s.ptype = obj.ptype;
            s.scale = obj.scale;
            s.problemID = obj.problemID;
            s.dirichlet = obj.dirichlet;
            s.pointload = obj.pointload;
            if isequal(obj.scale,'MICRO')
                s.masterSlave = obj.masterSlave;
            end
        end
        
    end
    
    methods (Access = private)
        
        function m = createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            s.masterSlaveNodes = obj.masterSlave;
            s.boundaryNodes    = obj.boundaryNodes;
            s.boundaryElements = obj.boundaryElements;
            m = Mesh(s);       
        end
        
        function readFile(obj,fileName)
            data = Preprocess.readFromGiD(fileName);
            [~,~,bNodes,bElem,mSlave] = Preprocess.getBC_mechanics(fileName);

            obj.pdim = data.problem_dim;
            obj.geometryType = data.geometry;
            obj.ptype = data.problem_type;
            obj.scale = data.scale;
            obj.problemID = fileName;
            obj.coord = data.xpoints;
            ndim = obj.getDimension();
            obj.coord  = obj.coord(:,2:ndim+1);
            obj.connec = data.connectivities(:,2:end);
            obj.boundaryNodes = bNodes;
            obj.boundaryElements = bElem;
            obj.masterSlave = mSlave;
            
            if strcmpi(data.problem_type,'elastic') ...
            || strcmpi(data.problem_type,'hyperelastic') ...
            || strcmpi(data.problem_type,'thermal')
                if isfield(data,'dirichlet_data')
                    obj.dirichlet = data.dirichlet_data;
                    obj.pointload = data.pointload;
                end
            end
        end
        
        function d = getDimension(obj)
            switch obj.pdim
                case '1D'
                    d = 1;
                case '2D'
                    d = 2;
                case '3D'
                    d = 3;
            end
        end
        
    end
    
end