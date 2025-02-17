classdef Preprocess<handle
    % Only reads
    
    %% !! SUSCEPTIBLE OF BEING RENAMED !!
    
    properties
    end
    
    methods(Static)
        
        function data = readFromGiD(filename)
            run(filename)
            data = struct;
            if exist('gidcoord','var')
                coord=gidcoord;
                connec=gidlnods;
            end
            data.xpoints=coord;
            if length(data.xpoints(1,:))==3
                data.xpoints(:,4)=0;
            end
            if any(connec(:,length(connec(1,:))))==0
                data.connectivities=connec(:,1:(length(connec(1,:))-1));
            else
                data.connectivities=connec;
            end
            data.geometry = strjoin(Data_prb(1));
            data.problem_dim = strjoin(Data_prb(3));
            data.problem_type = strjoin(Data_prb(5));
            data.scale = strjoin(Data_prb(6));
            
            if strcmpi(data.problem_type,'elastic')...
               || strcmpi(data.problem_type,'hyperelastic')...
               || strcmpi(data.problem_type,'thermal')
                if exist('dirichlet_data','var')
                    data.dirichlet_data = dirichlet_data;
                elseif exist('lnodes','var')
                    data.dirichlet_data = lnodes;
                else
                end
                
                if exist('pointload_complete','var')
                    data.pointload = pointload_complete;
                else
                    data.pointload = [];
                end
            end
        end
        
        function [fixnodes,forces,boundaryNodes,boundaryElements,Master_slave, sDir, sPL, sPer] = getBC_mechanics(filename)
            run(filename)
            if exist('lnodes','var')
                dirichlet_data=lnodes;
            end
            fixnodes{1,1} = dirichlet_data;
            
            if exist('pointload_complete','var')
                forces = pointload_complete;
            else
                forces = [];
            end
            
            if exist('External_border_nodes','var')
                boundaryNodes= External_border_nodes;
            else
                boundaryNodes = [];
            end

            if exist('External_border_elements','var')
                boundaryElements = External_border_elements;
            else
                boundaryElements = [];
            end            
            
            if ~isempty(boundaryNodes)
                boundaryNodes(:,2)=ones(length(boundaryNodes(:,1)),1);
                boundaryNodes(:,3)=zeros(length(boundaryNodes(:,1)),1);
            end
            
            if ~exist('Master_slave','var')
                Master_slave = [];
            end

            if exist('sDir','var')
                 sDir = sDir;
            else
                 sDir = [];
            end

            if exist('sPL','var')
                 sPL = sPL;
            else
                 sPL = [];
            end

            if exist('sPer','var')
                 sPer = sPer;
            else
                 sPer = [];
            end
            
        end
        
        function [state, velocity, pressure, Vol_force, velocityBC, dtime, finalTime, sDir] = getBCFluids(fileName)
            run(fileName)
            if ~exist('dtime', 'var')
                % Steady
                dtime = Inf;
                finalTime = [];
            end
            
        end
        
        function [forces_adjoint,pl] = getBC_adjoint(filename, mesh)
            run(filename)
            pl = [];
            for i = 1:numel(sPLAdj)
                pl = [pl, PointLoad(mesh, sPLAdj{i})];
            end
            forces_adjoint = pointload_adjoint;
        end
    end

end