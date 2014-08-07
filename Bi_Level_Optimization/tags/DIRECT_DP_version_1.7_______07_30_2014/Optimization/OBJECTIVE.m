function obj = OBJECTIVE(x,varargin)
obj=0;

% % Update the Design Variables
% assignin('base','FD',x(1));  % Why should I do this?

% Get_Params;
FD = x(1);
G = x(2);
fc_trq_scale = x(3);
mc_trq_scale = x(4);

% Run the model in the base workspace
Dynamic_Programming;

if FAIL == 1; % It failed
    FAIL_DP = 100;
else
    FAIL_DP = 0.5;
end

obj = -MPG; % UPDATE WITH EMMISSIONS/PERFORMANCE

assignin('base','con',FAIL_DP);

return
