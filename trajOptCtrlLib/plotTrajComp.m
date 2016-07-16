function h = plotTrajComp(t, trajectories, rows, cols, indicies, lineProp, plotTitle, labels, leg)
    % plotTraj(t, trajectories, rows, cols, indicies, lineProp, plotTitle, labels, leg)
    %   Simplifies comparison plots of trajectories using subplot()
    h = figure;
    for j = 1:rows*cols
        subplot(rows, cols, j);
        hold on; grid on; axis tight;
        title(labels{j});
        ylabel(labels{j});
        xlabel('t (s)');
        for k = 1:length(trajectories)
            % Should swap the nesting of these loops so we're not
            % constantly tranposing the same matrices
            [m, n] = size(trajectories{k});
            if m > n
                trajectories{k} = trajectories{k}';   % Make it wide
            end
            if(iscell(t))
                x = t{k};
            else
                x = t;
            end
            y = trajectories{k}(indicies(j), :);
            % Make sure vectors are the same length, truncate otherwise
            if length(x) > length(y)
                x = x(1:length(y));
            else
                y = y(1:length(x));
            end
            plot(x, y, lineProp{k});
        end
    end
    legend(leg);
end