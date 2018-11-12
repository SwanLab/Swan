function [trid] = find_point_element(conectivities,coordinatesa,gp_newmesh)
% A esta función le doy una pareja de puntos xy y me dice en que triangulo
% se encuentra.
% PASOS:
% 1) MALLA DE REFERENCIA conectivities,coordinatesa
% 2) IC; incentros de la malla de referencia
% 3) Cedges; lados de los triangulos
% 4) TRIA; malla de referencia; conectivities,coordinatesa
% 5) TRIA_HULL; triangulacion de los nodos de referencia con todos los lados
%    como restriccion. El resultado es el convexo con las restricciones.
% 6) etarget; elemento de TRIA_HULL que contiene incentro de TRIA  

TRIA = TriRep(conectivities', coordinatesa(:,1),coordinatesa(:,2));
IC = incenters(TRIA);
Cedges = edges(TRIA);
%triplot(TRIA);
% Triangulo el convexo
TRIA_HULL = DelaunayTri(coordinatesa, Cedges);
%seemesh('TRIA_HULL',1,coordinatesa,TRIA_HULL);
%triplot(TRIA_HULL);

% etarget; elementos de TRIA_HULL que contienen incentros de TRIA
etarget = pointLocation(TRIA_HULL,IC);
% lista de elementos de tria_hull; inicializada en nan
%dtTotrIndex = nan(size(TRIA_HULL,1),1);
dtTotrIndex = zeros(size(TRIA_HULL,1),1);
% dado elemento de tria_hull devuelve numero de elem de referencia o nan
dtTotrIndex(etarget) = 1:size(TRIA,1);

% trid; elemento de TRIA_HULL que contiene pto de gauss de malla nueva
trid = pointLocation(TRIA_HULL,gp_newmesh);
% elemento de TRIA que contiene pto de gauss de malla nueva 
if isfinite(trid)
    trid = dtTotrIndex(trid);
end
%keyboard
for i=1:size(trid,1)
    if trid(i)==0
       x=gp_newmesh(i,:);
       mindist=1e20;
       for j=1:size(IC,1)
%            dist = norm(x-IC(j)); 
           dist = norm(x-IC(j,:)); % corrección J.C. may 2013
           if dist<mindist
               mindist=dist;
               emin=j;
           end
       end
       trid(i)=emin;
    end
end



