clear all

S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.0002;
graph = 1;
[ V_0_new, V_f_new, Acc_Final_new ] = Accel_Req_2( S, dt_2, graph );  % V_Final is in MPH

