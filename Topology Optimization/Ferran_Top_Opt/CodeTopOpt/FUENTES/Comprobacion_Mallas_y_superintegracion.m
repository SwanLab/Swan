function Comprobacion_Mallas_y_superintegracion
paint = 0;
Vreal = compute_initial_radius(paint);
file_name = 'RVE04N3';
gaus = [1 3 4 6 7 13];
for imesh = 1:7
for igaus = 1:length(gaus)%nsnapshots
igaus
imesh
[~,~,~,~,~,~,~,~,~,vol] = read_data_problem(file_name,imesh,gaus(igaus)); 
volum(igaus,imesh) = vol;
error(igaus,imesh) = abs(vol - Vreal)/Vreal;
end
figure(2)
plot(gaus,error)
end

end

