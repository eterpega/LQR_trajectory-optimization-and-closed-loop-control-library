x0 = [30*pi/180 0]';
tspan = [0 50]
p.l1 = 1;
p.m1 = 2;
p.g = 9.81;
u = @(t) 1*sin(.5*2*pi*t);
[t, y] = ode45(@(t, x) pointMassPendulumEom(t, x, u(t), p), tspan, x0);

figure
plot(t, y(:, 1))