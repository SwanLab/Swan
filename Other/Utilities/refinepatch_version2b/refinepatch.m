function [FV2]=refinepatch(FV)
% This function "refinepatch" refines a triangular mesh with 
% a spline interpolated 4-split method.
%
%   [FV2] = refinepatch(FV,options)
%
% inputs,
%   FV : Structure containing a Patch, with
%        FV.vertices the mesh vertices
%        FV.face the mesh faces (triangles), rows with each 3 vertex indices
% outputs,
%   FV2 : Structure Containing the refined patch
%
%
% Reference:
%  The spline interpolation of the face edges is done by the 
%  Opposite Edge Method, described in: "Construction of Smooth Curves 
%  and Surfaces from Polyhedral Models" by Leon A. Shirman 
%
% How it works:
%  The tangents (normals) and velocity on the edge points of all edges 
%  are calculated. Which are  later used for b-spline interpolation when 
%  splitting the edges.
%
%  A tangent on an 3D line or edge is under defined and can rotate along 
%  the line, thus an (virtual) opposite vertex is used to fix the tangent and
%  make it more like a surface normal.
%
%  B-spline interpolate a half way vertices between all existing vertices
%  using the velocity and tangent from the edgepoints. After splitting a
%  new facelist is constructed
%
% Speed:
%  Compile the c-functions for more speed with:
%   mex vertex_neighbours_double.c -v;
%   mex edge_tangents_double.c -v;
%
% Example:
%
% X=[-0.5000;  0.5000;  0.0000;  0.0000];
% Y=[-0.2887; -0.2887;  0.5774;  0.0000];
% Z=[ 0.0000;  0.0000;  0.0000;  0.8165];
% FV.vertices=[X Y Z];
%
% FV.faces=[2 3 4; 4 3 1; 1 2 4; 3 2 1];
%
% figure, set(gcf, 'Renderer', 'opengl'); axis equal;
% for i=1:4
%   patch(FV,'facecolor',[1 0 0]);
%   pause(2);
%   [FV]=refinepatch(FV);
% end
%
% Function is written by D.Kroon University of Twente (February 2010)

% Get the neighbour vertices of each vertice from the face list.
Ne=vertex_neighbours(FV);

% Calculate the tangents (normals) and velocity of all edges. Which is
% later used for b-spline interpolation and split of the edges
%
% A tangent on an 3D line or edge is under defined and can rotate along 
% the line, thus an (virtual) opposite vertex is used to fix the tangent and
% make it more like a surface normal.
V=FV.vertices; F=FV.faces;
[ET_table,EV_table,ETV_index]=edge_tangents(V,Ne);

% B-spline interpolate a half way vertices between all existing vertices
% using the velocity and tangent from above
[V,HT_index, HT_values]=make_halfway_vertices(EV_table,ET_table,ETV_index,V,Ne);

% Make new facelist
Fnew=makenewfacelist(F,HT_index,HT_values);

FV2.vertices=V;
FV2.faces=Fnew;



