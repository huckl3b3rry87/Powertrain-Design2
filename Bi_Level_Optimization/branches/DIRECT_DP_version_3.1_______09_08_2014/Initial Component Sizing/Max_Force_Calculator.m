function [Sim_Grade, Force] = Max_Force_Calculator( vinf, dvar, V, Motor_ON  )

Wheel_Spd = V/vinf.rwh;   % rad/sec

x3_length = length(vinf.gear);
 
for x3 = 1:x3_length
    Gear = vinf.gear(x3);
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

