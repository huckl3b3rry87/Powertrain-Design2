function [Sim_Grade, Force] = Force_Calculator_Grade( vinf, dvar, V, Motor_ON  )

Wheel_Spd = V/vinf.rwh;   % rad/sec

x3_grid = [4.484 2.872 1.842 1.414 1.000 0.742];  % [1st 2nd...]             % Gear Level
x3_length = length(x3_grid);

for x3 = 1:x3_length
    Gear = x3_grid(x3);
    Eng_Spd = Wheel_Spd*dvar.FD*Gear;
    Eng_Torque =  interp1(vinf.eng_consum_spd,vinf.eng_max_trq,Eng_Spd);
    Eng_Wheel_Torque = Eng_Torque*dvar.FD*Gear;             %N*m
    Eng_Tractive_Effort(:,x3) = Eng_Wheel_Torque/vinf.rwh;  % N
end

Max_Eng_Tractive_Effort = max(Eng_Tractive_Effort,[],2); % Check the dim

% Motor
Motor_Spd = Wheel_Spd*dvar.FD*dvar.G;
Motor_Torque = interp1(vinf.m_map_spd,vinf.m_max_trq,Motor_Spd);
Motor_Wheel_Torque = Motor_Torque*dvar.FD*dvar.G;
Motor_Tractive_Effort = Motor_Wheel_Torque/vinf.rwh;

if Motor_ON == 1  % Include Motor Force
    Force = (Max_Eng_Tractive_Effort + Motor_Tractive_Effort);
    
else               % Engine Only
    Force = Max_Eng_Tractive_Effort;
end

Sim_Grade = [Eng_Tractive_Effort, Max_Eng_Tractive_Effort, Motor_Tractive_Effort, Force];

end

