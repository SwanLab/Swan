function printArrayInMatlabFormat(array)
    fprintf('[\n');
    for i = 1:size(array, 1)
        fprintf('    ',i);
        for j = 1:size(array, 2)
            fprintf('%i', array(i, j));
            if j < size(array, 2)
                fprintf(', ');
            end
        end
        if i < size(array, 1)
            fprintf(';\n');
        else
            fprintf('\n');
        end
    end
    fprintf('];\n');
end
