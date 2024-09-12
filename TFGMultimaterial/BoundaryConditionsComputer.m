classdef BoundaryConditionsComputer < handle

    properties (Access = public)
        bc
    end

    properties (Access = private)
        p
        e
        t
        mesh
        coord
        connec
    end

    methods (Access = public)

        function obj = BoundaryConditionsComputer(cParams)
            obj.init(cParams);
            obj.computeNeumannBoundaryConditions();
            obj.computeDirichletBoundaryConditions();
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.499*yMax & abs(coor(:,2))<=0.501*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = [1,2];
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.bc.pCons = []; %---> for static condensation
            obj.bc.pNeu = [];  
            obj.bc.Neu.ldir = [2.0,0.499,2.0,0.501];
            obj.bc.Neu.ldof = [1,2];
            obj.bc.Neu.lval = -1;
            obj.bc.pDir = [];
            obj.bc.Dir.ldir = [0.0,0.0,0.0,1.0];
            obj.bc.Dir.ldof = [2,1,2];
            obj.bc.Dir.lval = [0.0,0.0];

            obj.p = cParams.mesh.p;
            obj.e = cParams.mesh.e;
            obj.t = cParams.mesh.t; 

            obj.mesh = cParams.m; 

        end

        function computeNeumannBoundaryConditions(obj)
            s.ldir = obj.bc.Neu.ldir;
            s.ldof = obj.bc.Neu.ldof;
            s.lval = obj.bc.Neu.lval;

            [pbc] = computeBounCon(obj,s);
            obj.bc.pNeu = pbc;         
        end

        function computeDirichletBoundaryConditions(obj)
            s.ldir = obj.bc.Dir.ldir;
            s.ldof = obj.bc.Dir.ldof;
            s.lval = obj.bc.Dir.lval;
            
            pbc = computeBounCon(obj,s);
            obj.bc.pDir = pbc;
        end

        function pbc = computeBounCon(obj,s)
            pbc = []; % Initialize pbc here
            
            ldir = s.ldir;
            ldof = s.ldof;
            lval = s.lval;
            
            pb = obj.p;
            eb = obj.e;
            tb = obj.t;

            if ~isempty(ldir) && ~isempty(ldof) && ~isempty(lval)
                nlines = size(ldir,1);
                for j=1:nlines
                    [ndir,~,~] = findnodes(pb,eb,tb,ldir(j,:));
                    
                    if size(ndir,2)>1
                        ndir = ndir';
                    end
                    
                    ncc = size(ndir,1);
            
                    node = [];
                    dof = [];
                    val = [];
                    
                    for i=1:ldof(j,1)
                        node = [node;ndir];
                        dof  = [dof; ldof(j,i+1)*ones(ncc,1)]; 
                        val = [val;lval(j,i)*ones(ncc,1)];
                    end
            
                    bci = [node, dof, val];
                    pbc = [pbc;bci];
            
                end

            end

        end

    end

end