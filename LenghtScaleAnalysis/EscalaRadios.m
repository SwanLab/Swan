close all
clear all
mat = zeros(80,160);
poscionY = ceil(80/2);
poscionX = ceil(160/13);
radio = [1 2 3 4 5 6 7 8 9];
for rad = 1:length(radio)
    for elemX  = (poscionX-(radio(rad)+2)):(poscionX+(radio(rad)+2))
        for elemY  = (poscionY-(radio(rad)+2)):(poscionY+(radio(rad)+2))
            radElem = sqrt((elemX-poscionX)^2 + (elemY-poscionY)^2);
            if radElem < radio(rad)
                mat(elemY,elemX) = 1;
            end 
        end 
    end
    if rad ~= length(radio)
     poscionX = poscionX+radio(rad+1)*3;
    end
end 
mi_paleta = [1, 1, 1; 1, 0, 0];
figure(2)
colormap(mi_paleta)
imagesc(mat);

rutaGraf = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\','RadiousScale','.png'];
saveas(rutaGraf);
