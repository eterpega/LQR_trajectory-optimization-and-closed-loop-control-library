function doublePendCartTest
    close all;
%     x0 = [0; 150*pi/180; -150*pi/180; -1; 45*pi/180; 20*pi/180];
    x0 = [0; 150*pi/180; -150*pi/180; -.1; 0; 0];
%     p.m1 = 1; p.b1 = .1; 
%     p.m2 = .25; p.I2 = 5e-5; p.c2 = .3; p.l2 = .4; p.b2 = .01; 
%     p.m3 = .25; p.I3 = 5e-5; p.c3 = .2; p.l3 = .4; p.b3 = .01; 
%     p.g = 9.81;

    p.m1 = 0.125; p.b1 = 0;%0.01;
    p.m2 = 0.098; p.I2 = 0.00205; p.c2 = 0.20956; p.l2 = 0.422; p.b2 = 0.00125;
    % 2nd pendulum with bob at end
    % p.m3 = 0.13971; p.I3 = 0.0030865; p.c3 = 0.26557; p.l3 = 0.4; p.b3 = 0.005;
    % 2nd pendulum w/o bob
    p.m3 = 0.08646; p.I3 = 0.001446; p.c3 = 0.18; p.l3 = 0.411; p.b3 = 0.00125;
    p.g = 9.81;

    u = @(t, x) 0;%5*exp(-t)*sin(1*2*pi*t);
    [t, y] = ode45(@(t, x)doublePendCartEom(t, x, u(t, x), p), [0 40], x0);
    figure; hold on;
%     plot(t, y(:, 1), t, unwrap(y(:, 2)), t, unwrap(y(:, 2)));
    plot(t, y(:, 1), 'r');
    plot(t, y(:, 2), 'g');
    plot(t, y(:, 3), 'b');
    
    q1 = y(:, 1);
    q2 = y(:, 2);
    q3 = y(:, 3);
    legend('q_1', 'q_2', 'q_2');
%     return

    figure;
%     plot([-1 1], [0 0], 'k');
    xlim([-3 3]); ylim([-2 2]); axis equal; hold on;
    plot([-3 3], [0 0], 'k');
%     rod1Pos = @(q) [q(1), 0];
    % Anonymous function for position of origin of 2nd pendulum
    cartHeight = .25;
    cartWidth = .5;
    rod2Pos = @(x) [x(1)+p.l2*sin(x(2)), -p.l2*cos(x(2))];
    cart = createBox(-cartWidth/2, 0, .5, cartHeight, 0, [.9 .5 .1]);
    pend1 = createRod(q1(1), cartHeight/2, .035, p.l2, q2(1)-pi, 6, [.2 .2 .7]);
    rod2xy = rod2Pos([q1(1) q2(1) q3(1)]);
    pend2 = createRod(rod2xy(1), rod2xy(2) + cartHeight/2, .035, p.l3, q3(1)-pi, 6, [.9 .9 0]);
    tStart = tic;
    n=1;
    while n < length(t)
        tFrame = toc(tStart);           % Find elapsed time
%         [~, n] = min(abs(t-tFrame));    % Get corresponding frame
        while n < length(t) && t(n) < tFrame, n = n+1; end

        
%         tic
        rod2xy = rod2Pos([q1(n) q2(n) q3(n)]);
        movePoly(cart, q1(n)-cartWidth/2, 0, 0);
        movePoly(pend1, q1(n), cartHeight/2, q2(n)-pi);
        movePoly(pend2, rod2xy(1), rod2xy(2)+cartHeight/2, q3(n)-pi);
%         toc
        xlim([-3 3]);
        ylim([-2 2]);
        title(sprintf('t: %2.1fs', t(n)));
%         drawnow;
        pause(.01);
    end
end

function p = createBox(x, y, width, height, angle, colour)
    p.geometry = [[0 0 width width]', [0 height height 0]'];
    p.vertices = rotate([p.geometry(:,1)+x p.geometry(:,2)+y], angle);   
    p.col = colour;
    p.x = x; p.y = y; p.angle = angle;
    p.h = fill(p.vertices(:, 1), p.vertices(:, 2), colour, 'EdgeColor', 'none');
end

function rotv = rotate(vertices, theta)
    [m, n] = size(vertices);
    for k = 1:m
        rotv(k, :) = ([cos(theta) -sin(theta); sin(theta) cos(theta)]*vertices(k, :)')';
    end
end

% Create a polygon and display it
function p = createRod(x, y, width, height, angle, nCapPoints, colour)
    cap = width*[cos(linspace(0, pi, nCapPoints))', sin(linspace(0, pi, nCapPoints))'];
    p.geometry = [cap(:, 1) cap(:,2)+height; -cap(:, 1) -cap(:, 2)];
    p.vertices = rotate([p.geometry(:,1)+x p.geometry(:,2)+y], angle);
    p.col = colour;
    p.x = x; p.y = y; p.angle = angle;
    p.h = fill(p.vertices(:, 1), p.vertices(:, 2), colour, 'EdgeColor', 'none');
end

function movePoly(poly, x, y, angle)
    poly.vertices = rotate(poly.geometry, angle);
    set(poly.h, 'Vertices', [poly.vertices(:,1)+x poly.vertices(:,2)+y], 'LineSmoothing', 'on');
end    