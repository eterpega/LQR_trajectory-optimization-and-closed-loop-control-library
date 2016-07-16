function p = createRod(x, y, x0, y0, width, height, angle, nCapPoints, colour)
    % createRod(x, y, x0, y0, width, height, angle, nCapPoints, colour) 
    % Creates a polygon representation of a rod with round ends.
    % By default the rod is vertically orientated, with the origin at the
    % base - x0 and y0 shift the origin relative to this point. All further
    % rotations and translations will be relative to this origin.
    %
    %   x, y           position of origin
    %   x0, y0         displace origin relative to base (described above)
    %   width, height  Size of rod
    %   nCapPoints     Number of points for curve
    %   color          Letter or vector ('r', [1 0 0] etc)
    %   p              Handle to created polygon object
    
    cap = width*[cos(linspace(0, pi, nCapPoints))', sin(linspace(0, pi, nCapPoints))'];
    p.geometry = [cap(:,1)-x0 cap(:,2)+height-y0; -cap(:, 1)-x0 -cap(:, 2)-y0];
    p.vertices = rotatePoly([p.geometry(:,1) p.geometry(:,2)], angle);
    p.vertices = [p.vertices(:, 1)+x, p.vertices(:, 2)+y];
    p.col = colour;
    p.x = x; p.y = y; p.angle = angle;
    p.h = fill(p.vertices(:, 1), p.vertices(:, 2), colour, 'EdgeColor', 'none');
end