classdef SplitterInConnectedComponents < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        faces
    end
    
    methods (Access = public)
        
        function obj = SplitterInConnectedComponents(cParams)
            obj.init(cParams)            
        end
        
        function C = split(obj)
            C = obj.connected_components(obj.faces);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.faces = cParams.faces;
        end
        
        function C = connected_components(obj,F)
            A = obj.adjacency_matrix(F);
            [~,C] = obj.conncomp(A);
        end
        
        function [A] = adjacency_matrix(obj,E)
            if size(E,2)>2
                F = E;
                E = obj.meshEdges(F);
            end
            A = sparse([E(:,1) E(:,2)],[E(:,2) E(:,1)],1);
        end
        
        function [S,C] = conncomp(obj,G)
            % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
            % toobox. G is an n by n adjacency matrix, then this identifies the S
            % connected components C. This is also an order of magnitude faster.
            %
            % [S,C] = conncomp(G)
            %
            % Inputs:
            %   G  n by n adjacency matrix
            % Outputs:
            %   S  scalar number of connected components
            %   C
            
            % Transpose to match graphconncomp
            G = G';
            
            [p,~,r] = dmperm(G+speye(size(G)));
            S = numel(r)-1;
            C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
            C(p) = C;
        end
        
        function edges = meshEdges(obj,faces, varargin)
            if isstruct(faces) && isfield(faces, 'faces')
                % if input is a mesh structure, extract the 'faces' field
                faces = faces.faces;
            elseif nargin > 2
                % if two arguments are given, keep the second one
                faces = varargin{1};
            end
            
            
            if ~iscell(faces)
                % Process faces given as numeric array
                % all faces have same number of vertices, stored in nVF variable
                
                % compute total number of edges
                nFaces  = size(faces, 1);
                nVF     = size(faces, 2);
                nEdges  = nFaces * nVF;
                
                % create all edges (with double ones)
                edges = zeros(nEdges, 2);
                for i = 1:nFaces
                    f = faces(i, :);
                    edges(((i-1)*nVF+1):i*nVF, :) = [f' f([2:end 1])'];
                end
                
            else
                % faces are given as a cell array
                % faces may have different number of vertices
                
                % number of faces
                nFaces  = length(faces);
                
                % compute the number of edges
                nEdges = 0;
                for i = nFaces
                    nEdges = nEdges + length(faces{i});
                end
                
                % allocate memory
                edges = zeros(nEdges, 2);
                ind = 0;
                
                % fillup edge array
                for i = 1:nFaces
                    % get vertex indices, ensuring horizontal array
                    f = faces{i}(:)';
                    nVF = length(f);
                    edges(ind+1:ind+nVF, :) = [f' f([2:end 1])'];
                    ind = ind + nVF;
                end
                
            end
            
            % keep only unique edges, and return sorted result
            edges = sortrows(unique(sort(edges, 2), 'rows'));
            
        end
        
        
        
    end
    
    
end