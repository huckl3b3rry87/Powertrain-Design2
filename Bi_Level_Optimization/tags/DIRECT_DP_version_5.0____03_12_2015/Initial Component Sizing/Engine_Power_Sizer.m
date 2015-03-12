function [Sim_Grade, F_max_6, V_max_6, RR ] = Engine_Power_Sizer( param, vinf, dvar, test )

% Desired Grade and Top Speed Performance - (Engine) - For Low Grades
V = [0:0.1:120*param.mph_mps]';

%----------------------Set Requirements ---------------------------------
% Test 1
r = 1;
V_test(r) = 32.5*param.mph_mps;
alpha_test(r) = 0*pi/180;

% Test 2
r = 2;
V_test(r) = 20*param.mph_mps;
alpha_test(r) = 6*pi/180;

% -------------------- Base Tests ---------------------------------
RR = length(alpha_test);
for r = 1:1:RR
    % Pbase(r) = Engine_Base_Power(V_test(r),alpha_test(r), param, vinf);
    % P_req(r,:) = Engine_Base_Power(V,alpha_test(r), param, vinf);
    F_RL(:,r) = Road_Load( V, alpha_test(r), param, vinf ); 
end

% ---------------- Actual Engine Performance -----------------------------
Motor_ON = 0;

[Sim_Grade, Force_PM] = Max_Force_Calculator( vinf, dvar, V, Motor_ON  );

% For all gears
[F_max_t, V_max_t ] = Max_Speed( Force_PM, F_RL(:,test), V );

w = length(vinf.gear);
[ F_max_6, V_max_6] = Max_Speed(Sim_Grade(:,w), F_RL(:,test), V );

%Sim_Grade = [Eng_Tractive_Effort, Max_Eng_Tractive_Effort, Motor_Tractive_Effort, Force_PM];
Sim_Grade = [Sim_Grade, F_RL, V, V_max_t*ones(size(V)), F_max_t*ones(size(V))];

end

