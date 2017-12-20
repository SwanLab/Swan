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
        
        function [fixnodes, forces] = getBC(filename)
            run(filename)
            fixnodes = lnodes;
            forces = pointload_complete;
        end
        
        function [fixnodes, pnods] = getPeriodicBC(coordinates)
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
            
            nodes = 1:length(coordinates(:,1));
            %Vertex 1 y 2
            index_bin_max_y = coordinates(:,2) == max(coordinates(:,2));
            index_max_y = nodes(index_bin_max_y);
            
            [~,index1_local] = min(coordinates(index_max_y,1));
            index_vertex_1 = index_max_y(index1_local);
            [~,index2_local] = max(coordinates(index_max_y,1));
            index_vertex_2 = index_max_y(index2_local);
            
            %Vertex 3 y 4
            index_bin_min_y = coordinates(:,2) == min(coordinates(:,2));
            index_min_y = nodes(index_bin_min_y);
            
            [~,index3_local] = max(coordinates(index_min_y,1));
            index_vertex_3 = index_min_y(index3_local);
            [~,index4_local] = min(coordinates(index_min_y,1));
            index_vertex_4 = index_min_y(index4_local);
            
            index_vertex = [index_vertex_1;index_vertex_2;index_vertex_3;index_vertex_4];
            corners = index_vertex';
            
            % prescription of the corners
            ifix=0;
            for i=1:size(corners,2)
                ifix=ifix+1;
                fixnodes(ifix,1)=corners(i); % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=corners(i); % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            end
        end
    end
end





