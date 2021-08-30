classdef FEM < handle
    
    properties (GetAccess = public, SetAccess = public)
        problemData
        geometry
        mesh
        dof
        element
        Fext
        variables
    end
    
    properties (Access = protected)
        solver
        iter = 0;   
        inputReader
    end
    
    methods (Access = public)
        
        function obj = FEM()
            obj.inputReader = FemInputReader_GiD();
        end
        
    end
    
    methods (Static, Access = public)
        
        function obj = create(fileName)
            s = FemInputReader_GiD().read(fileName);
            switch s.ptype
                case 'ELASTIC'
                    switch s.scale
                        case 'MACRO'
                            obj = Elastic_Problem(fileName);
                        case 'MICRO'
                            obj = Elastic_Problem_Micro(fileName);
                    end
                case 'THERMAL'
                    obj = Thermal_Problem(fileName);
                case 'DIFF-REACT'
                    obj = DiffReact_Problem(fileName);
                case 'HYPERELASTIC'
                    obj = Hyperelastic_Problem(fileName);
                case 'Stokes'
                    obj = Stokes_Problem(fileName);
            end
        end
        
    end
    
    methods (Access = public)
        
        function sol = solve_steady_nonlinear_problem(obj,free_dof,tol)
            total_free_dof = sum(free_dof);
            dr = obj.element.computedr;
            x0 = zeros(total_free_dof,1);
            
            r = obj.element.computeResidual(x0,dr);
            x = x0;
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r
                r = obj.element.computeResidual(x,dr);
                x0 = x;
            end
            sol = x;
        end
        
        function sol = solve_transient_nonlinear_problem(obj,free_dof,tol,dt,final_time)
            total_free_dof = sum(free_dof);
            x_n(:,1) = zeros(total_free_dof,1);
            x0 = zeros(total_free_dof,1);
            
            dr = obj.element.computedr(dt);
            
            for istep = 2: final_time/dt
                u_previous_step = x_n(1:free_dof(1),istep-1);
                
                r = obj.element.computeResidual(x0,dr,u_previous_step);
                while dot(r,r) > tol
                    inc_x = obj.solver.solve(dr,-r);
                    x = x0 + inc_x;
                    % Compute r
                    r = obj.element.computeResidual(x,dr,u_previous_step);
                    x0 = x;
                end
                x_n(:,istep) = x;
            end
            sol = x_n;
        end
        
        function setDof(obj,dof)
            obj.dof = dof;
        end
        
        function setMatProps(obj,s)
           obj.element.material.compute(s);           
        end
        
        function setC(obj,C)
            obj.element.material.C = C;
        end
        
        function print(obj,fileName)
            dI = obj.createPostProcessDataBase(fileName);
            postprocess = Postprocess('Elasticity',dI);
            q = obj.element.quadrature; 
            d.fields = obj.variables;
            d.quad = q;
            postprocess.print(obj.iter,d);
        end
        
        function print_slave(obj,~,evnt)
            postprocess = Postprocess_PhysicalProblem;
            res_file = evnt.AffectedObject.res_file;
            postprocess.print_slave(obj,res_file,obj.variables);
        end        
        
       function syncPostProcess(obj,evtobj)
            addlistener(evtobj,'res_file','PostSet',@obj.print_slave);
        end        
        
        function i = getIter(obj)
            i = obj.iter;
        end
        
    end
    
    methods (Access = protected)
        
        function readProblemData(obj,fileName)
            femReader = FemInputReader_GiD();
            s = femReader.read(fileName);
            
            obj.problemData.fileName = fileName;
            obj.problemData.scale = s.scale;
            obj.problemData.pdim  = s.pdim;
            obj.problemData.ptype = s.ptype;
            obj.problemData.nelem = s.mesh.nelem;
            obj.problemData.bc.dirichlet = s.dirichlet;
            obj.problemData.bc.pointload = s.pointload;
            obj.mesh = s.mesh;
        end
        
        function createMesh(obj)
            s.coord  = obj.inputReader.coord;
            s.connec = obj.inputReader.connec;            
            obj.mesh = Mesh(s);
        end
        
    end
    
    methods (Access = private)
        
        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh    = obj.mesh;
            dI.outName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.getValue();           
        end
        
    end
    
    methods (Access = public, Abstract)
        %preProcess(obj)
        computeVariables(obj)
    end
    
end
