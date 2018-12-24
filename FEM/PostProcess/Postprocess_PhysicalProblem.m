classdef Postprocess_PhysicalProblem < Postprocess
    

    
    methods (Access = public)
      
        
        function  print(obj,physical_problem,file_name,physicalVars)
            results.physicalVars = physicalVars;
            path = pwd;
            dir = fullfile(path,'Output',file_name);
            if ~exist(dir,'dir')
                mkdir(dir)
            end
            obj.setBasicParams(physical_problem,file_name,results)
            obj.PrintMeshFile();
            obj.PrintResFile(results)
        end
        
        
    end
    
    methods (Access = protected)
        function setBasicParams(obj,physical_problem,file_name,results)
            obj.nfields = physical_problem.element.nfields;
            for ifield = 1:obj.nfields
                obj.coordinates{ifield} = physical_problem.element.interpolation_u(ifield).xpoints;
                obj.connectivities{ifield} = physical_problem.element.interpolation_u(ifield).T;
                obj.nnode(ifield) = physical_problem.element(ifield).nnode;
                obj.npnod(ifield) = physical_problem.element.interpolation_u(ifield).npnod;  % Number of nodes
            end
            obj.gtype = physical_problem.mesh.geometryType;
            obj.ndim = physical_problem.element.interpolation_u.ndime;
            obj.pdim = physical_problem.mesh.pdim;
            obj.ngaus = physical_problem.element(1).quadrature.ngaus;
            obj.posgp = physical_problem.element(1).quadrature.posgp';
            obj.ptype = physical_problem.mesh.ptype;
            
            switch  obj.gtype % GiD type
                case 'TRIANGLE'
                    obj.etype = 'Triangle';
                case 'QUAD'
                    obj.etype = 'Quadrilateral';
                case 'TETRAHEDRA'
                    obj.etype = 'Tetrahedra';
                case 'HEXAHEDRA'
                    obj.etype = 'Hexahedra';
            end
            obj.nelem = physical_problem.element.nelem; % Number of elements
            
            obj.gauss_points_name = 'Guass up?';
            
            obj.file_name = file_name;
            switch obj.ptype
                case 'ELASTIC'
                    obj.nsteps = 1;
                case 'Stokes'
                    obj.nsteps = length(results.physicalVars.u(1,:));
            end
        end
    end
    
    methods (Access = private)
        
        function PrintResFile(obj,results)
            iD = obj.fid_res;
            fN = obj.file_name;
            nS = obj.nsteps;
            gD = obj.gauss_points_name;
            eT = obj.etype;
            pT = obj.ptype;
            nG = obj.ngaus;
            nD = obj.ndim;
            pG = obj.posgp;
            rS = results;
            ElasticityResultsPrinter(iD,fN,nS,gD,eT,pT,nG,nD,pG,rS);
        end
        
        
    end
end

