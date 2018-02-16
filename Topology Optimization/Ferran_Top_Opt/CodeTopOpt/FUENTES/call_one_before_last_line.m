function lastLine2 = call_one_before_last_line(fid)

lastLine = '';                   %# Initialize to empty
offset = 1;                      %# Offset from the end of file
fseek(fid,-offset,'eof'); %# Seek to the file end, minus the offset
newChar = fread(fid,1,'*char');  %# Read one character
while (~strcmp(newChar,char(10))) || (offset == 1)
    lastLine = [newChar lastLine];   %# Add the character to a string
    offset = offset+1;
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    newChar = fread(fid,1,'*char');  %# Read one character
end

lastLine2 = '';
newChar = fread(fid,1,'*char');
while (~strcmp(newChar,char(10))) || (offset == 1)
    lastLine2 = [newChar lastLine2];   %# Add the character to a string
    offset = offset+1;
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    newChar = fread(fid,1,'*char');  %# Read one character
end

end