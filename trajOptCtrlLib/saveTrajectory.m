function saveTrajectory(matFile)
% loadTrajectory(matFile) loads the trajectory from a previous trajectory
% optimization saved in matFile and resamples it to fit nGridPoints grid
% points. 
% Horrible way to do this, but need to think about how I'm going to package
% optimized trajectory solutions before I sort this out.
fStr = ['save(''' matFile ''', ''traj'', ''u'', ''T'', ''param'', ''output'');'];
evalin('caller', fStr);