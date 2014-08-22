function [Max_Motor_Spd ] = Max_Motor_Spd( Motor_Tractive_Effort, V )

A =  isnan(Motor_Tractive_Effort);  % If there is a one then there is a NaN
I = find(A==0);  % if a never = 0, all nan..
if isempty(I)
    pause;
end
Max_Motor_Spd = V(I(end));

end

