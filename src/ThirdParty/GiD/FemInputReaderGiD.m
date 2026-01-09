classdef FemInputReaderGiD < handle
    
    properties (GetAccess = public,SetAccess = public)
        pdim
        dirichlet
        pointload
        mesh

        dirichletFun
        pointloadFun
        periodicFun
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
            s.mesh = obj.mesh;
            s.pdim = obj.pdim;
            s.geometryType = obj.geometryType;
            s.ptype = obj.ptype;
            s.scale = obj.scale;
            s.problemID = obj.problemID;
            s.dirichlet = obj.dirichlet;
            s.pointload = obj.pointload;
            s.dirichletFun = obj.dirichletFun;
            s.pointloadFun = obj.pointloadFun;
            s.periodicFun  = obj.periodicFun;
            if isequal(obj.scale,'MICRO')
                s.masterSlave = obj.masterSlave;
            end
            if isequal(obj.ptype,'Stokes')
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
            m = Mesh.create(s);
        end
        
        function readFile(obj,fileName)
            data = Preprocess.readFromGiD(fileName);
            if isequal(data.problem_type,'Stokes')
                [state, vel, prs, forces, velBC, dtime, fTime,sDir] = Preprocess.getBCFluids(fileName);
                obj.state = state;
                obj.velocity = vel;
                obj.pressure = prs;
                obj.velocityBC = velBC;
                obj.forcesFormula = forces;
                obj.dtime = dtime;
                obj.finalTime = fTime;
                % sDir = [];
                sPL = [];
                sPer = [];
            end
            if ~isequal(data.problem_type,'Stokes')
                [~,~,bNodes,bElem,mSlave, sDir, sPL, sPer] = Preprocess.getBC_mechanics(fileName);
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
            obj.mesh = obj.createMesh();
            
            obj.dirichletFun = [];
            obj.pointloadFun = [];
            obj.periodicFun  = [];
            if ~isequal(sDir,[])
                for i = 1:numel(sDir)
                    dir = DirichletCondition(obj.mesh, sDir{i},data.dirichlet_data);
                    obj.dirichletFun = [obj.dirichletFun, dir];
                end
            end

            if ~isequal(sPL,[])
                for i = 1:numel(sPL)
                    pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                    obj.pointloadFun = [obj.pointloadFun, pl];
                end
            end

            if ~isequal(sPer,[])
                for i = 1:numel(sPer)
                    per = PeriodicCondition(obj.mesh, sPer{i});
                    obj.periodicFun = [obj.periodicFun, per];
                end
            end

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