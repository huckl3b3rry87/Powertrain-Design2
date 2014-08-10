% Program to determine Fixed Throttle Performance (FTP)

% make this program easy to use and simulate everything - maybe a gui..
% different engines and motors
% What is your Acceleration time? - That is what this program can do
% What is your top speed? - That is what this program can do
% What is your desired grade?  - Cannot do this yet...
% What is your engine, motor, vehicle..ect.
% Scale the engine up based off first equation
% Determine size of motor based off this
% Then identify a battery
% Outputs are all of the scaling factors and the number of modules
% Later you can do a better job selecting gear ratios

clc; clear all; close all;
global m rw ni g grade Cd Af rho Frr eng_max_trq eng_consum_spd W_eng_min W_eng_max
% Share information with the ‘function’
mph_mps = 1/2.237;
% Enter proper vehicle parameters:
m= 1.9013e+03; rw=0.287; Frr=0.015;
g = 9.81; grade = 0; Cd = 0.28; Af = 2.52; rho = 1.2; 
eng_max_trq= 3*[61 67.6 73.7 78.5 80.9 77.3 76.2 73.3 68.7]; 
eng_consum_spd=[104.5 149.2 220.9 292.5 364.1 435.7 507.4 552.2 596.9];
W_eng_min = min(eng_consum_spd);
W_eng_max = max(eng_consum_spd); 

% Gearbox and final drive ratios:
ng= [4.484 2.872 1.842 1.414 1.000 0.742]; nf=3.5;
x0=[eps 0]; t0=0; tf = 0.5; % Initial conditions
wem=eng_consum_spd(end-2); % Engine speed at which gear is shifted
we=0; % A value set for the ‘while’ statement below
for i=1: length(ng) % Loop for gears
    ni=ng(i)*nf; % Overall gear ratio at each gear
    while we < wem % Repeat statements until we=wem
        [t,x]=ode45(@Fixed_thrt, [t0 tf], x0); % Calls Fixed_thrt function
        v=x(:,1);  % Speed
        s=x(:,2);  % Distance
              
        We_c = ni*v/rw; % rad/sec
        We_c(We_c < W_eng_min) = W_eng_min;
        We_c(We_c > W_eng_max) = W_eng_max;
        Trq = interp1(eng_consum_spd,eng_max_trq,We_c);
        Fti=Trq*ni/rw;
        Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*v.^2*Af;
        acc=(Fti-Frl)/m;% Differential equation for speed
        if mean(acc) < 0.05
            break;
        end
        
        we = ni*v(end)/rw;
        if we < wem
            tf = tf + 0.5;
        end
    end
    we_sim = 30*ni*v/rw/pi; % RPM
    figure(3)
    plot(t, we_sim);
    hold on
    
    figure(4);
    plot(t,Trq);
    hold on
    
    % At this point results for gear ‘ni’ are ready to plot
    figure(1);
    subplot(2,1,1), plot(t,v/mph_mps), hold on
    subplot(2,1,2), plot(t,s), hold on
    % Now gearshift is needed. Set the initial conditions for the next gear
    t0=max(t);   

    tf = tf + 0.5; 
    x0=[v(end) s(end)];
    we=0;
    
    We_c = ni*v/rw; % rad/sec
    We_c(We_c < W_eng_min) = W_eng_min;  
    We_c(We_c > W_eng_max) = W_eng_max;
    Trq = interp1(eng_consum_spd,eng_max_trq,We_c);
    Fti=Trq*ni/rw;
    Frl = m*g*sin(grade) + Frr*m*g*cos(grade) + 0.5*rho*Cd*(v).^2*Af;
    
    acc=(Fti-Frl)/m;% Differential equation for speed
    figure(2)
    plot(t, acc')
    hold on
    clear acc
end
figure(1)
subplot(2,1,1)
xlabel('Time (s)')
ylabel('Velocity MPH')
grid
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Distance (m)')
grid

figure(2)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
grid

figure(3)
xlabel('Time (s)')
ylabel('Engine Speed (RPM)')
grid

figure(4)
xlabel('Time (s)')
ylabel('Engine Torque (Nm)')
grid

