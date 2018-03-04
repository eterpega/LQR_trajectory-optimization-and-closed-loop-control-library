% traj = x_traj_cl;
% t = t0;

traj = x_traj_cl;
t = linspace(0, T, nPoints)

[m, n] = size(traj);

% t = linspace(0, T, n);
f = figure;
hold on; axis equal; grid on;
xlimits = [-3 3];
xlim(xlimits);
ylim([-1 1.5]); 

plot(xlimits, [0 0], 'k');
cart = createBox(0, 0, .15, 0, .3, .25, 0, [.9 .5 .1]);
pend1 = createRod(traj(1, 1), .125, 0, 0, .025, sys.param.l2, traj(2, 1), 12, [.9 .9 0]);
pend2 = createRod(traj(1, 1), .125+sys.param.l2, 0, 0, .025, sys.param.l3, traj(3, 1), 12, [.1 .5 .9]);

trail.num = 0;
trail.delay = .3;
trail.col = [.7 .7 .7];
filename = '';%'dpendanim';
animTraj(traj', t, sys.param, { {cart, {'x', 'traj(n, 1)'}}, ...
                                {pend1, {'x', 'traj(n, 1)'}, {'ang', 'traj(n, 2)+pi'}}, ...
                                {pend2, {'x', 'traj(n, 1)+p.l2*sin(traj(n, 2))'}, ... 
                                        {'y', '.125-p.l2*cos(traj(n, 2))'}, ... 
                                        {'ang', 'traj(n, 3)+pi'}} ...
                }, trail, 1, filename);