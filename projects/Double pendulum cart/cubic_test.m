% cubic test function
% function cubic_test

clear vars;

% x = 0:.1:10;%[0 1 2.5 3.6 5 7 8.1 10];
% y = sin(x);
% xx = 0:.1:.5;
% yy = spline(x,y,xx);
% figure;
% plot(x,y,'o',xx,yy)

% % Some data
% t = 0.3:0.01:0.5;
% dt = t(end) - t(end-1);
% f = sin(2*pi*0.5*t);
% df = diff(f)/dt;

% Analytical version
syms fa dfa ta
fa = @(t) sin(t);
dfa = @(t) cos(t);
t = 0:.1:6;
ti = [];
fi = [];
steps = 10;
for ts = linspace(0, t(end)-t(end)/steps, steps)
    t0 = ts;
    t1 = ts+t(end)/steps;
    f0 = fa(t0);
    f1 = fa(t1);
    df0 = dfa(t0);
    df1 = dfa(t1);
    fi(end+1) = f0;
    ti(end+1) = t0;
    % Do some interpolating
    step = 10;
    dt = t1-t0;
    for k = 1/step:1/step:1
        ti(end+1) = t0+(k*dt);
        fi(end+1) =    (df1 - 2*f1 + df0 + 2*f0)*k^3 + ...
                (3*f1 - df1 - 2*df0 - 3*f0)*k^2 + df0*k + f0;
%         fi(end+1) =    (df1 - 2*f1 + df0 + 2*f0)*k^3 + ...
%                 (3*f1 - df1 - 3*f0)*k^2 + df0*k + f0;
    end
end
figure; hold on; grid on;
plot(t, subs(fa, t), '-.');
plot(ti, fi);
plot(linspace(0, t(end), steps+1), fa(linspace(0, t(end), steps+1)), 'xk')
return



% Some other data
t = 0:.1:3;
f = sin(t);%sin(2*pi*0.1*t);
dt = t(end) - t(end-1);
df = diff(f)/dt;
% df = [df(1) df];

% % Still more data
% 
% t = 0:1:10;
% f = sin(t);
% dt = t(end) - t(end-1);
% df = diff(f)/dt;


figure;
plot(t, f, '.', t, df, '-.');
grid on
% df(end+1) = df(end);

% fi, ti, interpolated t and f vectors
ti = [];
fi = [];
dfi =[];

ngrid = [1 2];%[5 6];%
% ngrid = 1:length(t);
dt = t(ngrid(2)) - t(ngrid(1));
for n=1:length(ngrid)-1%1:length(t)-1
    t0 = t(ngrid(n));
    f0 = f(ngrid(n));
    f1 = f(ngrid(n+1));
    df0 = df(ngrid(n));
    df1 = df(ngrid(n+1));

    fi(end+1) = f0;
    ti(end+1) = t0;
    % Do some interpolating
    step = 10;
    for k = 1/step:1/step:1
        ti(end+1) = t0+(k*dt);
        fi(end+1) =    (df1 - 2*f1 + df0 + 2*f0)*k^3 + ...
                (3*f1 - df1 - 2*df0 - 3*f0)*k^2 + df0*k + f0;
    end
    dfi = diff(fi)/((t(ngrid(2)) - t(ngrid(1)))/step);
end
hold on
plot(ti, fi)
plot(linspace(t(ngrid(1)), t(ngrid(2)), step), spline(t, f, linspace(t(ngrid(1)), t(ngrid(2)), step)))
% plot(linspace(t(ngrid(1)), t(ngrid(2)), step), dfi);