function PrintFileINFO(FILE_BATCH,TEXTP)

    fid =fopen(FILE_BATCH,'w');
    for i = 1:length(TEXTP)
        fprintf(fid,[TEXTP{i},'\n']);
    end
    fod =fclose(fid);
    
    open(FILE_BATCH);