function mesh = mesh_square_X_solid(L, hole_scale, Nr, Ntheta)
    % L: Mitad del lado del cuadrado (ej. 1 para un cuadrado de 2x2)
    % hole_scale: Tamaño del hueco (0.6 = barras gruesas, 0.9 = barras finas)
    % Nr: Divisiones entre el hueco y el borde
    % Ntheta: Divisiones a lo largo de cada lado del hueco
    
    % 1. Definir el centro y las 4 esquinas del cuadrado
    centro = [0, 0];
    esquinas = [ L,  L;   % 1. Superior Derecha
                -L,  L;   % 2. Superior Izquierda
                -L, -L;   % 3. Inferior Izquierda
                 L, -L];  % 4. Inferior Derecha
    
    nodes_all = [];
    elements_all = [];
    node_offset = 0;

    % 2. Generar 4 sectores (uno por cada triángulo que apunta al centro)
    for k = 1:4
        % Vértices del triángulo del sector (Externos)
        v1_ext = centro;
        v2_ext = esquinas(k,:);
        v3_ext = esquinas(mod(k,4)+1,:);
        
        sector_ext = [v1_ext; v2_ext; v3_ext];
        
        % Calcular baricentro del sector para escalar el hueco
        baricentro = mean(sector_ext);
        
        % Vértices del hueco (Internos)
        % Se escalan respecto al baricentro para dejar espacio en todos los bordes
        v1_int = baricentro + hole_scale * (v1_ext - baricentro);
        v2_int = baricentro + hole_scale * (v2_ext - baricentro);
        v3_int = baricentro + hole_scale * (v3_ext - baricentro);
        sector_int = [v1_int; v2_int; v3_int];

        % 3. Mallar el área sólida (3 trapezoides por sector)
        % Cada trapezoide conecta un lado del triángulo interno con el externo
        for lado = 1:3
            idx_next = mod(lado,3)+1;
            
            p_int1 = sector_int(lado,:);
            p_int2 = sector_int(idx_next,:);
            p_ext1 = sector_ext(lado,:);
            p_ext2 = sector_ext(idx_next,:);
            
            [n_trap, e_trap] = create_trapezoid_mesh(p_int1, p_int2, p_ext1, p_ext2, Nr, Ntheta);
            
            % Acumular nodos y elementos
            e_trap = e_trap + node_offset;
            nodes_all = [nodes_all; n_trap];
            elements_all = [elements_all; e_trap];
            node_offset = size(nodes_all, 1);
        end
    end

    % 4. Fusionar nodos duplicados (Crucial para que la X y el marco sean una sola pieza)
    tol = 1e-8;
    [nodes_final, ~, ic] = unique(round(nodes_all/tol)*tol, 'rows', 'stable');
    xmax = L;
    xmin = -L;
    ymax = L;
    ymin = -L;
    elements_final = ic(elements_all);

                mesh.coord = nodes_final;
            mesh.coord(mesh.coord(:,1)== xmax & mesh.coord(:,2)==ymax,:) =...
                mesh.coord(mesh.coord(:,1)== xmax & mesh.coord(:,2)==ymax,:)+[-1e-9,-1e-9];
            mesh.coord(mesh.coord(:,1)== xmax & mesh.coord(:,2)==ymin,:) =...
                mesh.coord(mesh.coord(:,1)== xmax & mesh.coord(:,2)==ymin,:)+[-1e-9,1e-9];
            mesh.coord(mesh.coord(:,1)== xmin & mesh.coord(:,2)==ymax,:) =...
                mesh.coord(mesh.coord(:,1)== xmin & mesh.coord(:,2)==ymax,:)+[1e-9,-1e-9];
            mesh.coord(mesh.coord(:,1)== xmin & mesh.coord(:,2)==ymin,:) =...
                mesh.coord(mesh.coord(:,1)== xmin & mesh.coord(:,2)==ymin,:)+[1e-9,1e-9];

    % Organizar salida
%     mesh.coord = nodes_final;
    mesh.connec = elements_final;

%     % Visualización
%     figure;
%     patch('Faces', mesh.connec, 'Vertices', mesh.coord, ...
%           'FaceColor', [0.8 0.9 1], 'EdgeColor', [0.2 0.2 0.2]);
%     axis equal; grid on;
%     title(['Marco sólido con X interna (hole\_scale = ', num2str(hole_scale), ')']);
end

function [nodes, elements] = create_trapezoid_mesh(int1, int2, ext1, ext2, Nr, Ntheta)
    % Crea una malla de cuadriláteros entre dos segmentos de recta
    nodes = [];
    for j = 0:Ntheta
        t = j / Ntheta;
        % Interpolación en el borde del hueco y en el borde exterior
        pA = (1 - t) * int1 + t * int2;
        pB = (1 - t) * ext1 + t * ext2;
        
        for i = 0:Nr
            s = i / Nr;
            % Interpolación radial (del interior al exterior)
            nodes = [nodes; (1 - s) * pA + s * pB];
        end
    end

    % Generar tabla de conectividades
    node_grid = reshape(1:size(nodes,1), Nr+1, Ntheta+1);
    elements = [];
    for j = 1:Ntheta
        for i = 1:Nr
            % Definir los 4 nodos del elemento cuadrangular
            elements(end+1, :) = [node_grid(i,j), node_grid(i+1,j), ...
                                  node_grid(i+1,j+1), node_grid(i,j+1)];
        end
    end
end