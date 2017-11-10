% Incremental encoder test program

syms t;
x = @(t) [5*sin(t*2*pi); cos(0.5*t*2*pi); 10*pi*cos(t*2*pi); -pi*sin(0.5*t*2*pi)];
controlFreq = 20;   % Hz
encoderLines = 1440;
% incEnc_old(0, x(0), controlFreq, encoderLines, false);
enc = IncEnc(2, 0, controlFreq, encoderLines, true);
tv = 0:0.01:2;
x_traj_c = x(tv);

x_traj_d = zeros(4, length(tv));
for n = 1:length(tv)
%     x_traj_d(:, n) = incEnc_old(tv(n), x(tv(n)));    
    x_traj_d(:, n) = enc.read(tv(n), x(tv(n)));
end
% figure
% plot(tv, x_traj_c(1, :), 'b', tv, x_traj_c(2, :), 'm');
% hold on
% plot(tv, x_traj_d(1, :), 'b--', tv, x_traj_d(2, :), 'm--');
% legend('Pos (''continuous'' time)', 'Speed (''continuous'' time)', 'Pos estimate', 'Speed estimate');
% xlabel('t');
% title(sprintf('Controller at %d Hz', controlFreq));

figure; hold on;
plot(enc.Trajectory.xIn.t, enc.Trajectory.xIn.traj(1,:), 'b', enc.Trajectory.xIn.t, enc.Trajectory.xIn.traj(2,:), 'c');
plot(enc.Trajectory.xIn.t, enc.Trajectory.xIn.traj(3,:), 'g', enc.Trajectory.xIn.t, enc.Trajectory.xIn.traj(4,:), 'y');

plot(enc.Trajectory.xEst.t, enc.Trajectory.xEst.traj(1,:), 'b--', enc.Trajectory.xEst.t, enc.Trajectory.xEst.traj(2,:), 'c--');
plot(enc.Trajectory.xEst.t, enc.Trajectory.xEst.traj(3,:), 'g--', enc.Trajectory.xEst.t, enc.Trajectory.xEst.traj(4,:), 'y--');
% legend('Pos (''continuous'' time)', 'Speed (''continuous'' time)', 'Pos estimate', 'Speed estimate');
xlabel('t');
title(sprintf('Controller at %d Hz', controlFreq));