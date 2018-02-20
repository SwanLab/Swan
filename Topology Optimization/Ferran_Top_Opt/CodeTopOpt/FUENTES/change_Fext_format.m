
function New_Fext = change_Fext_format(Fext)

nnode = length(unique(Fext(:,1)));
New_Fext = zeros(nnode,3);

for inode = 1:nnode
    New_Fext(inode,1) = Fext(inode,1);
    node = Fext(:,1) == Fext(inode,1);
    New_Fext(inode,Fext(node,2)+1) = Fext(node,3);
end

end