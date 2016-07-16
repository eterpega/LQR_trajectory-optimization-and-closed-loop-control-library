function p = createBox(x, y, x0, y0, width, height, angle, colour)
    p.geometry = [[0 0 width width]'-x0, [0 height height 0]'-y0];
    % Need to fix this so it works like createRod:
    p.vertices = rotatePoly([p.geometry(:,1)+x p.geometry(:,2)+y], angle);   
    p.col = colour;
    p.x = x; p.y = y; p.angle = angle;
    p.h = fill(p.vertices(:, 1), p.vertices(:, 2), colour, 'EdgeColor', 'none');
end