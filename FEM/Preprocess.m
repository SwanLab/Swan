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
        
        function [fixnodes,fixnodes_perimeter, forces] = getBC(filename)
            run(filename)
            fixnodes = lnodes;
            forces = pointload_complete;
            fixnodes_perimeter=External_border_nodes;
            if ~isempty(fixnodes_perimeter)
                fixnodes_perimeter(:,2)=ones(length(fixnodes_perimeter(:,1)),1);
                fixnodes_perimeter(:,3)=zeros(length(fixnodes_perimeter(:,1)),1);
            end
        end
        function forces_adjoint=getBC_adjoint(filename)
            run(filename)
            forces_adjoint = pointload_adjoint;
        end
        
        function [pnods] = getPeriodicBC(coordinates)
            % PERIODIC BOUNDARY COND
            % creation of a list containing the couples that define the periodicity
            % 1) list of nodes on each side (two vertical, two horizontal)
            % 2) sort each list based on the corresponding coordinate
            % 3) generation of the couples and store them in 'pnods'
            %
            % remark:
            % [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
            % If X is a vector, then Y = X(I).
            
            % nodes in the left-vertical side, without the corners
            href = 0.025;
            h=href; L=[];Y=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,1)<h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5 )
                    L = [L i];
                    Y = [Y coordinates(i,2)];
                end
            end
            [Y1,I] = sort(Y);
            V1 = L(I);
            
            % nodes in the right-vertical side, without the corners
            h=href; L=[];Y=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,1)>1-h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5)
                    L = [L i];
                    Y = [Y coordinates(i,2)];
                end
            end
            [Y1,I] = sort(Y);
            V2 = L(I);
            
            % nodes in the bottom-horizontal side, without the corners
            h=href; L=[];X=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,2)<h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5 )
                    L = [L i];
                    X = [X coordinates(i,1)];
                end
            end
            [X1,I] = sort(X);
            H1 = L(I);
            
            % nodes in the top-horizontal side, without the corners
            h=href; L=[];X=[];
            for i=1:size(coordinates,1)
                if (coordinates(i,2)>1-h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5)
                    L = [L i];
                    X = [X coordinates(i,1)];
                end
            end
            [X1,I] = sort(X);
            H2 = L(I);
            pnods = [V1 H1; V2 H2]; % lista de nodos
        end
    end
end





