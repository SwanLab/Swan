classdef Preprocess<handle
    % Only reads
    
    %% !! SUSCEPTIBLE OF BEING RENAMED !!
    
    properties
    end
    
    methods(Static)
        function data = readFromGiD(filename)
            
            addpath(fullfile('.','Input'))
            run(filename)
            data = struct;
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
        end
        
        function [fixnodes,forces,full_dirichlet_data,Master_slave] = getBC(filename)
            run(filename)
            fixnodes = dirichlet_data;
            forces = pointload_complete;
            
            full_dirichlet_data=External_border_nodes;
            if ~isempty(full_dirichlet_data)
                full_dirichlet_data(:,2)=ones(length(full_dirichlet_data(:,1)),1);
                full_dirichlet_data(:,3)=zeros(length(full_dirichlet_data(:,1)),1);
            end
            
            if ~exist('Master_slave','var')
               Master_slave = []; 
            end
            
        end
        function forces_adjoint=getBC_adjoint(filename)
            run(filename)
            forces_adjoint = pointload_adjoint;
        end
        

    end
end





