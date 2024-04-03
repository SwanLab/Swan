classdef plotejarmalla

    methods(Static)
        function plotPunts(punts) %s.mesh.coord
            % Obtenció de les coordenades x i y de la matriu
            x = punts(:, 1);
            y = punts(:, 2);
            
            % Plot dels punts amb unitats entre ells
            plot(x, y, 'o');
            xlabel('x');
            ylabel('y');
            title('');

            for i = 1:size(punts, 1)
                text(x(i), y(i), num2str(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            end

        end
    end



end
