clear variables;
% derive_dpendcart();
load('doublePendCartSys.mat', 'sys');

sys.name = 'doublePendCartTest';

% If parameters change, run these two lines to update gradient function
sys = createDircolNlConGrads2(sys);
save([sys.name 'Sys.mat'], 'sys');
