% encProperties - vector of encoder resolutions (post quad, so true
% resolution encoder will see. Zero to bypass
function xOut = incEnc_old(t, x, freqInit, encRes, recordTraj)
    persistent freq;        % Controller frequency
    persistent tLast;
    persistent encLast;
    persistent encLines;
    persistent nStates;
    persistent xEst;
    persistent recordOn;
    persistent history;
    
    if(nargin > 2)
        % Initialise
        freq = freqInit;
        tLast = t;
        encLines = encRes;
        nStates = length(x)/2;
        encLast = zeros(nStates, 1);
        xEst = zeros(length(x), 1);
        if(recordTraj == true)
            % Record the trajectory
            recordOn = true;
            history.xin.t = [];
            history.xin.traj = [];
            history.xest.t = [];
            history.xest.traj = [];
        end
        return
    end
    
    if(encLines == 0)
        xOut = x;
        return;
    else
        if t >= tLast+1/freq
            encNow = floor((x(1:nStates)/(2*pi)).*encLines);
            % Time to update
            enc_delta = encNow - encLast;
            enc_dt = (enc_delta./encLines)*2*pi/(t-tLast);
            xEst = [2*pi*encNow./encLines; enc_dt];
            encLast = encNow;
            tLast = t;
%             if(recordOn)
                
        end
        xOut = xEst;
    end
%     if(recordOn)
%         history.xin.t = [history.xin.t t];
%         history.xin.traj = [history.xin.traj x];
%     end
end
