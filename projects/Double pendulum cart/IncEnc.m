classdef IncEnc < handle
    properties
        Freq;
        TLast;
        EncLast;
        EncLines;
        NumVars;
        XEst;
        Record;
        Trajectory;
    end
    
    methods
        function obj = IncEnc(nVars, tStart, frequency, encoderResolution, recordTrajOnOff)
            obj.Freq = frequency;
            obj.TLast = tStart;
            obj.EncLines = encoderResolution;
            obj.NumVars = nVars;
            obj.EncLast = zeros(nVars, 1);
            obj.XEst = zeros(nVars*2, 1);
            if(recordTrajOnOff)
                % Record the trajectory
                obj.Record = true;
                obj.Trajectory.xIn.t = [];
                obj.Trajectory.xIn.traj = [];
                obj.Trajectory.xEst.t = [];
                obj.Trajectory.xEst.traj = [];
            end
        end
        
        function xOut = read(obj, t, x)
            % Bypass?
            if(obj.EncLines == 0)
                xOut = x;
            else
                if t >= obj.TLast+1/obj.Freq
                    % Time to update
                    encNow = floor((x(1:obj.NumVars)/(2*pi)).*obj.EncLines);
                    enc_delta = encNow - obj.EncLast;
                    enc_dt = (enc_delta./obj.EncLines)*2*pi/(t-obj.TLast);
                    obj.XEst = [2*pi*encNow./obj.EncLines; enc_dt];
                    obj.EncLast = encNow;
                    obj.TLast = t;
%                     if(obj.Record)
%                         obj.Trajectory.xEst.t = [obj.Trajectory.xEst.t t];
%                         obj.Trajectory.xEst.traj = [obj.Trajectory.xEst.traj obj.XEst];
%                     end
                end
                xOut = obj.XEst;
            end
                        % Store input trajectory
            if(obj.Record)
                obj.Trajectory.xIn.t = [obj.Trajectory.xIn.t t];
                obj.Trajectory.xIn.traj = [obj.Trajectory.xIn.traj x];
                obj.Trajectory.xEst.t = [obj.Trajectory.xEst.t t];
                obj.Trajectory.xEst.traj = [obj.Trajectory.xEst.traj obj.XEst];
            end
        end
        
        function getTrajectory(obj)
            % Actually not needed because we can access class properties
            % directly
            
            
            
            
        end
    end
end