% Given matrix
matrix = [
    1 1 1;
    1 1 2;
    1 1 3;
    1 2 1;
    1 2 2;
    1 3 1;
    2 1 1;
    2 1 2;
    2 2 1;
    3 1 1
];

% Define a mapping function
mapping = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:size(matrix, 1)
    key = mat2str(matrix(i, :));
    if ~isKey(mapping, key)
        mapping(key) = i;
    end
end

% Test the mapping function
input1 = [1 3 1];
input2 = [3 1 1];

output1 = mapping(mat2str(input1));
output2 = mapping(mat2str(input2));

disp(['Mapping for ' mat2str(input1) ' is: ' num2str(output1)]);
disp(['Mapping for ' mat2str(input2) ' is: ' num2str(output2)]);
