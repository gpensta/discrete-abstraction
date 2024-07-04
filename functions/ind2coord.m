function [i, j] = ind2coord(g, n)
    % Convert unique index g to (i, j) coordinates
    % n is the number of columns in the grid
    i = floor((g - 1) / n) + 1;
    j = g - (i - 1) * n;
end