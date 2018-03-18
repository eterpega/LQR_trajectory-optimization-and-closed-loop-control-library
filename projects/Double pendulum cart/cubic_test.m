% Cubic Hermite spline (I believe)
% https://dsp.stackexchange.com/questions/18265/bicubic-interpolation/18273#18273
% cubic test function
% function cubic_test

% clear vars;

% Another go (continuous functions)
f = @(t) sin(t);
df = @(t) cos(t);
t = 0:.5:10;

ti = [];
fi = [];
dfi = [];
steps = 10;
for ts = linspace(0, t(end)-t(end)/steps, steps)
    t0 = ts;
    t1 = ts+t(end)/steps;
    f0 = f(t0);
    f1 = f(t1);
    df0 = df(t0);
    df1 = df(t1);
    fi(end+1) = f0;
    ti(end+1) = t0;
    dfi(end+1) = df0;
    % Do some interpolating
    kstep = 5;
    dt = t1-t0;
    for k = 1/kstep:1/kstep:1
        k
        t0+(k*dt)
        ti(end+1) = t0+(k*dt);
        fi(end+1) =    (df1 - 2*f1 + df0 + 2*f0)*k^3 + ...
                (3*f1 - df1 - 2*df0 - 3*f0)*k^2 + df0*k + f0;
        dfi(end+1) = 3*(df1 - 2*f1 + df0 + 2*f0)*k^2 + 2*(3*f1 - df1 - 3*f0 - 2*df0)*k + df0;
%         fi(end+1) =    (df1 - 2*f1 + df0 + 2*f0)*k^3 + ...
%                 (3*f1 - df1 - 3*f0)*k^2 + df0*k + f0;
    end
end
figure; hold on; grid on;
plot(t, subs(f, t), '-.');
plot(ti, fi);
% plot(ti, dfi);
% plot(t, df(t));
plot(linspace(0, t(end), steps+1), f(linspace(0, t(end), steps+1)), 'xk')
return

% Some other data
% t = 1:length(lqr.K(:,1));
f = lqr.K(:,1)';
dt = t(end) - t(end-1);
df = diff(f)/dt;
df = [df df(end)];

figure;
plot(t, f, '.');%, t, df, '-.');
grid on
% df(end+1) = df(end);

% fi, ti, interpolated t and f vectors
ti = [];
fi = [];
% dfi =[];

% ngrid = [1 2];%[5 6];%
ngrid = 1:20:length(t);
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
%     dfi = diff(fi)/((t(ngrid(2)) - t(ngrid(1)))/step);
end
hold on
plot(ti, fi)
% plot(linspace(t(ngrid(1)), t(ngrid(2)), step), spline(t, f, linspace(t(ngrid(1)), t(ngrid(2)), step)))
% plot(linspace(t(ngrid(1)), t(ngrid(2)), step), dfi);