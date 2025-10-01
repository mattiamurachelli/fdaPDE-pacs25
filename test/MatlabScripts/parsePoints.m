% Function to convert a column of "x1;x2;..." strings into a numeric matrix
function mat = parsePoints(PointCol)
    n = height(PointCol);          % number of rows
    % Preallocate based on first row
    firstRow = split(PointCol{1}, ';');
    d = length(firstRow);             % dimension of the vector
    mat = zeros(n, d);
    
    for i = 1:n
        parts = split(PointCol{i}, ';');   % split string by ';'
        mat(i,:) = str2double(parts)';        % convert to numeric row vector
    end
end
