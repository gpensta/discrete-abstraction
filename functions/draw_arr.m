function draw_arr(x, y)
    headSize = .05;
    dx = y(1) - x(1);
    dy = y(2) - x(2);
    mag = sqrt(dx^2 + dy^2);
    dx_unit = dx / mag;
    dy_unit = dy / mag;
    dx_short = dx - headSize * dx_unit;
    dy_short = dy - headSize * dy_unit;
    h1 = quiver(x(1), x(2), dx_short, dy_short, 0, 'k', 'LineWidth', .5, 'MaxHeadSize', 0);
    hold on;
    h2 = quiver(x(1) + dx_short, x(2) + dy_short, headSize * dx_unit, headSize * dy_unit, 0, 'k', 'LineWidth', .5, 'MaxHeadSize', 2);
    if nargout > 0
        varargout{1} = [h1, h2];
    end
end