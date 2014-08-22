clear V_0 V_f
S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.0002;
graph = 1;
[ V_0, V_f, Acc_Final ] = Accel_Req_2( S, dt_2, graph );  % V_Final is in MPH

for i = 1:length(Acc_Final)
[ PASS, Sim_Variables ] = Acceleration_Test_2(V_0(i),V_f(i), Acc_Final(i), dt_2, param, vinf, dvar);
PASS_sim(i)= PASS
end

%%
% i = 2;
% V_0(i)
% V_f(i)