classdef FemInputReader_GiD < handle
    
    properties (GetAccess = public,SetAccess = public)
        pdim
        dirichlet
        pointload
        ptype
        scale
        problemID
        fileName

        velocity
        pressure
        forcesFormula
        velocityBC
        state
        dtime
        finalTime
        nu
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
            if isequal(obj.ptype,'Stokes')
                s.nu       = obj.nu;
                s.state    = obj.state;
                s.dtime    = obj.dtime;
                s.ftime    = obj.finalTime;
                s.velocity = obj.velocity;
                s.pressure = obj.pressure;
                s.forcesFormula = obj.forcesFormula;
                s.velocityBC    = obj.velocityBC;
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
            if isequal(data.problem_type,'Stokes')
                [preData] = Preprocess.getBCFluidsNew(fileName);
                obj.nu = preData.nu;
                obj.state = preData.state;
                obj.velocity = preData.velocity;
                obj.pressure = preData.pressure;
                obj.velocityBC = preData.velocityBC;
                obj.forcesFormula = preData.Vol_force;
                obj.dtime = preData.dtime;
                obj.finalTime = preData.finalTime;
            end
            if ~isequal(data.problem_type,'Stokes')
                [~,~,bNodes,bElem,mSlave] = Preprocess.getBC_mechanics(fileName);
                obj.boundaryNodes = bNodes;
                obj.boundaryElements = bElem;
                obj.masterSlave = mSlave;
            end

            obj.pdim = data.problem_dim;
            obj.geometryType = data.geometry;
            obj.ptype = data.problem_type;
            obj.scale = data.scale;
            obj.problemID = fileName;
            obj.coord = data.xpoints;
            ndim = obj.getDimension();
            obj.coord  = obj.coord(:,2:ndim+1);
            obj.connec = data.connectivities(:,2:end);

            
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