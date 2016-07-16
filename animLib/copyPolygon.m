function p = copyPolygon(old, col)
    % p = copyPolygon(old, col)
    % Creates a new copy of a polygon, identical except for colour (col)
    p.h = fill(old.vertices(:, 1), old.vertices(:, 2), col, 'EdgeColor', 'none');
    p.geometry = old.geometry;
    p.vertices = old.vertices;
    p.col = col;
    p.x = old.x; p.y = old.y; p.angle = old.angle;
end