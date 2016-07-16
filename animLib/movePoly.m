function movePoly(poly, x, y, angle)
%       Translate and rotate a polygon
%           poly - handle to polygon
%           x - x coordinate
%           y - y coordinate
%           angle - angle of rotation (radians)
    poly.vertices = rotatePoly(poly.geometry, angle);
    set(poly.h, 'Vertices', [poly.vertices(:,1)+x poly.vertices(:,2)+y], 'LineSmoothing', 'on');
end    