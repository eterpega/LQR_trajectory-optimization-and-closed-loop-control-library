function animTraj(traj, t, p, coordVars, trail, speed, filename)
    % animTraj(traj, t, p, coordVars)
    %   traj        trajectory matrix (wide)
    %   t           time vector
    %   p           system physical parameters
    %   coordVars   expressions for calculating polygon coordinates

    if nargin > 5
        t = t/speed;
    end
    if nargin > 4 && isstruct(trail)
        trailson = true;
        delay = trail.delay;
        numTrails = trail.num;
        col = trail.col;
    else 
        trailson = false;
    end
    
    if trailson
        trails = {};
        for tNum = 1:numTrails
            for n = 1:length(coordVars)
                trails{end+1} = copyPolygon(coordVars{n}{1}, col); %#ok<AGROW>
                uistack(trails{end}.h, 'bottom');
            end
            col = col+([1 1 1]-col)*(tNum/(numTrails+1));
        end
    end

    exportVideo = nargin > 6 && ~strcmp(filename, '');
    if exportVideo
        vid = VideoWriter(filename, 'MPEG-4');
        % Doesn't support variable frame rate, so grab dt in middle of
        % animation and use to calculate frame rate.
        skip = 1; % If too many frames
        n = round(length(t)/4);
        vid.FrameRate = 1/((t(n+1)-t(n))*skip);
        open(vid);
    end
    
    tStart = tic;
    n=0;
    fr = 1;
    while n < length(t)
        if ~exportVideo % Don't save
            % Play in real time
            tFrame = toc(tStart);           % Find elapsed time
            [~, n] = min(abs(t-tFrame));    % Get corresponding frame
        else
            % Saving, so get every frame
            n = n + skip;
            if n > length(t), n=length(t); end
            tFrame = t(n);
        end
        for j = 1:length(coordVars)
            x=coordVars{j}{1}.x; y=coordVars{j}{1}.y; ang=coordVars{j}{1}.angle;
            for k=2:length(coordVars{j})
                eval([coordVars{j}{k}{1} '=' coordVars{j}{k}{2} ';']);
            end
            movePoly(coordVars{j}{1}, x, y, ang);
        end
%         title(sprintf('t: %2.2fs', t(n)));
        title(sprintf('t: %2.2fs', n));

        % Trails
        if trailson
            for tNum = 1:numTrails
                tTrail = tFrame - delay*tNum;
                if tTrail < 0
                    n = 1;
                else
                    [~, n] = min(abs(t-tTrail)); 
                end
                for j = 1:length(coordVars)
                    x=coordVars{j}{1}.x; y=coordVars{j}{1}.y; ang=coordVars{j}{1}.angle;
                    for k=2:length(coordVars{j})
                        eval([coordVars{j}{k}{1} '=' coordVars{j}{k}{2} ';']);
                    end
                    movePoly(trails{j+(tNum-1)*length(coordVars)}, x, y, ang);
                end
            end
        end
        % Update display
        drawnow;
        if exportVideo %nargin > 6 && ~strcmp(filename, '')
            % Based on animated GIF example from http://uk.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
            frame = getframe(gca);
            writeVideo(vid, frame);
%             im = frame2im(frame);
%             [imind, cm] = rgb2ind(im, 256);
%             if fr ~= length(t)
%                 fDuration = t(fr+1) - t(fr);
%             else
%                 fDuration = 3;
%             end
%             fDuration = 0.05*5;
%             if fr == 1
%                 imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
%             else
%                 imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append'); %, 'DelayTime', fDuration
%             end
%             fr
%             fr = fr + 1;
        end
    end
    if exportVideo, close(vid); end
end