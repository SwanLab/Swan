%% ********** CONVERSION FROM STRUCTURED TO UNSTRUCTURED MESH *************
function [A1,b1,A0,b0] = conversionTensors3D(problemID,dim,div)
run(problemID);

%     % Display
%     figure, hold on
%     for i = 1:50
%         if mod(i,1)==0 || i == 1
%             text(coord(i,2),coord(i,3),num2str(coord(i,1)));
%             plot(coord(i,2),coord(i,3),'b.')
%         end
%     end
%     axis('equal')

% Create index matrix and vector - P1
A1 = zeros(div(3)+1,div(2)+1,div(1)+1);
for n = 1:length(coord)
    inode = coord(n,1); x = coord(n,2); y = coord(n,3); z = coord(n,4);
    %         i = 1 + div(3) - round((z/dim(3))*div(3));
    %         j = 1 + div(2) - round((y/dim(2))*div(2));
    %     k = 1 + div(1) - round((x/dim(1))*div(1));
    
    i = 1 + round((z/dim(3))*div(3));
    j = 1 + round((y/dim(2))*div(2));
    k = 1 + round((x/dim(1))*div(1));
    
    A1(i,j,k) = inode;
    b1(inode,1) = i;
    b1(inode,2) = j;
    b1(inode,3) = k;
end

% Create index matrix and vector - P0
A0 = zeros(div(3),div(2),div(1));
for n = 1:length(connec)
    ielem = connec(n,1);
    for m = 1:4
        x_(m) = coord(connec(n,m+1),2);
        y_(m) = coord(connec(n,m+1),3);
        z_(m) = coord(connec(n,m+1),4);
    end
    x = mean(x_); y = mean(y_);  z = mean(z_);
    elem_coord(n,1) = x;
    elem_coord(n,2) = y;
    elem_coord(n,3) = z;
    
    %     i = div(3) - round((z/dim(3))*(div(3)-1));
    %     j = div(2) - round((y/dim(2))*(div(2)-1));
    %     k = div(1) - round((x/dim(1))*(div(1)-1));
    
    i = 1 + round((z/dim(3))*(div(3)-1));
    j = 1 + round((y/dim(2))*(div(2)-1));
    k = 1 + round((x/dim(1))*(div(1)-1));
    
    A0(i,j,k) = ielem;
    b0(ielem,1) = i;
    b0(ielem,2) = j;
    b0(ielem,3) = k;
end
end