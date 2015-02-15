function [Sim_Grade, FAIL_GRADE_TEST, V_max] = Grade_Test( param, vinf, dvar, alpha_test, V_test, Motor_ON )

% Desired Grade and Top Speed Performance - (Engine) - For Low Grades
V = [0:0.1:120*param.mph_mps]';

% -------------------- Base Tests ---------------------------------
RR = length(alpha_test);
for r = 1:1:RR
    F_RL(:,r) = Road_Load( V, alpha_test(r), param, vinf ); 
end

% ---------------- Actual Engine Performance -----------------------------

[Sim_Grade, Force_PM] = Max_Force_Calculator( vinf, dvar, V, Motor_ON );

if isempty(Force_PM)
    FAIL_GRADE_TEST = 1;
    V_max = [0 0];
else
    for r = 1:1:RR
        [F_max_t, V_max_t ] = Max_Speed( Force_PM, F_RL(:,r), V );
        
        if V_max_t < V_test(r)
            FAIL(r) = 1;
            V_max(r) = 0;
        else
            FAIL(r) = 0;
            V_max(r) = V_max_t;
                        
        end
        
    end
    
    %Sim_Grade = [Eng_Tractive_Effort, Max_Eng_Tractive_Effort, Motor_Tractive_Effort, Force_PM];
    Sim_Grade = [Sim_Grade, F_RL, V, V_max_t*ones(size(V)), F_max_t*ones(size(V))];
    
    FAIL_GRADE_TEST = any(FAIL);
end

end

