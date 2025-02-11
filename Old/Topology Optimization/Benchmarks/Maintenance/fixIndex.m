function new_name = fixIndex(name,destination_folderpath)
    index = 1;
    list = updateList(destination_folderpath);
    
    try_name = assembleName(name,index);
    
    for i = 1:length(list)
        if strcmpi(try_name,list(i).name)
            index = index + 1;
            try_name = assembleName(name,index);
        end
    end
    new_name = try_name;
end

function new_name = assembleName(name,index)
    us = strfind(name,'_');
    dots = strfind(name,'.');
    name(us(end)+1:dots(end)-1)=[];
    new_name = [name(1:us(end)), num2str(index), name(us(end)+1:end)];
end