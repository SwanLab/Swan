function move_inside=CheckMoveInside3D(V,Fnew,Vertex_id)

% Orignal vertex position
P=V(Vertex_id,:);
    
move_inside=false(1,size(Fnew,1));
for i=1:size(Fnew,1)
    A=V(Fnew(i,1),:); B=V(Fnew(i,2),:); C=V(Fnew(i,3),:);
    normal = -cross(A-C,B-C); 
    normal=normal./sqrt(sum(normal.^2));

    % Calculate shortest distance from vertice point to new face (plane)
    C=PointToClosestPointOnPlane(A,B,C,P);

    % Shortest distance vector
    normal2=P-C; normal2=normal2./sqrt(sum(normal2.^2));

    % Moves inside if they are opposite (equation equal to zero)
    move_inside(i)=sum((normal-normal2).^2)<1e-5;
end

% All face move inside, check ...
move_inside = min(move_inside);
