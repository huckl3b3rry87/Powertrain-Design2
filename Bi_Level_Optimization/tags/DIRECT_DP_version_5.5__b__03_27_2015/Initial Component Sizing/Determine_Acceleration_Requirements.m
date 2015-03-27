clear all
close all
S = 10;                    % Number of points to calcualte acceleration at
dt_2 = 0.0002;
graph = 1;
[ V_0, V_f, Acc_Final] = Accel_Req_Ann_Arbor( S, dt_2, graph );  % V_Final is in MPH

V_0 = V_0(1:7);
V_f = V_f(1:7);
Acc_Final = Acc_Final(1:7);

save('V_0_aa','V_0')
save('V_f_aa','V_f')
save('Acc_Final_aa','Acc_Final')