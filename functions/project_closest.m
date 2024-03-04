function [closest_x, closest_y, row, col] = project_closest(xgrid, ygrid, step_size, x, y)
    row = round((y - ygrid(1)) / step_size) + 1;
    col = round((x - xgrid(1)) / step_size) + 1;
    row = max(row, 1);
    col = max(col, 1);
    row = min(row, size(ygrid, 1));
    col = min(col, size(xgrid, 2));
    closest_x = xgrid(row, col);
    closest_y = ygrid(row, col);
end