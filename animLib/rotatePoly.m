function rotv = rotatePoly(vertices, theta)
    [m, n] = size(vertices);
    for k = 1:m
        rotv(k, :) = ([cos(theta) -sin(theta); sin(theta) cos(theta)]*vertices(k, :)')';
    end
end