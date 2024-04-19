function h = multimat_plot( p,t,fi )
%multimat_plot Plots the caracteristic function fi of each material and
%also the void region. Plots up to 4 different materials and void.

%auth: Augusto Romero
%date: 26-06-2018

n = size(fi,2);
colors = [0 0 0;1 0 0; 0 0 1;0 1 0];
h = gcf;

for i=1:n
    h(i) = pdesurf(p,t,fi(:,i));
    hold on,
    h(i).FaceColor = colors(i,:);
end
axis image; axis off;
view(2)
h(end).FaceColor = [1 1 1];
zlim([1e-1 1])

