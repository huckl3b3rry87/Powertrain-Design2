clear
clc

a2 = 0;
a1 = 0;

for k = 1:5
    
    if k ==1
        a3 = 0.1;
    elseif k ==2
        a3 = 0.25;
    elseif k == 3
        a3 = 1;
    elseif k == 4
        a3 =5;
    else
        a3 = 10;
    end
    
    [ FAIL, MPG, emission, delta_SOC ] = find_weights(a1,a2,a3)
    
    result.mpg(k) = MPG;
    result.NOx(k) = emission.NOx;
    result.HC(k)= emission.HC;
    result.CO(k) = emission.CO;
    result.dSOC(k) = delta_SOC;
    
    
end
Total = [result.mpg; result.NOx; result.HC; result.CO; result.dSOC]