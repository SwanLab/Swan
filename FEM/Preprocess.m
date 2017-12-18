classdef Preprocess<handle
    % Only reads
    
   %% !! SUSCEPTIBLE OF BEING RENAMED !! 
    
    properties
    end
    
    methods(Static)
        function data = readFromGiD(filename)
            addpath('.\Input\')
            run(filename)
            data = struct;
            data.xpoints=gidcoord;
            if length(data.xpoints(1,:))==3
                data.xpoints(:,4)=0;
            end
            if any(gidlnods(:,length(gidlnods(1,:))))==0
                data.connectivities=gidlnods(:,1:(length(gidlnods(1,:))-1));
            else
                data.connectivities=gidlnods;
            end
            data.geometry = strjoin(Data_prb(1));
            data.problem_dim = strjoin(Data_prb(3));
            data.problem_type = strjoin(Data_prb(5));
            data.scale = strjoin(Data_prb(6));
        end
        
        function [fixnodes, forces] = getBC(filename)
            run(filename)
            fixnodes = lnodes;
            forces = pointload_complete;
        end
        
    end
end



