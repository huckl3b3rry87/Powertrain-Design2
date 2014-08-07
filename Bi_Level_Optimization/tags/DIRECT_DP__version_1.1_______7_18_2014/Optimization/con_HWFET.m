function [con, con_e]=con_HWFET(x,varargin) 

con = evalin('base','con');

% offset=length(con);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Run Acceleration Test - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     cyc_mph = [0:100; 0:100]';
%     cyc_mph(:,2) = 200;
%     time_final = length(cyc_mph(:,2));
%     
%     V_prof=cyc_mph;
%     time_shift=1;
%     V_prof_pre=[cyc_mph(1+time_shift:end,1); zeros(time_shift,1)];
%     
%     init_soc = 0.65;
%     assignin('base', 'ess_init_soc', init_soc);
%     
%     assignin('base', 'V_prof',  V_prof);
%     assignin('base', 'time_shift', time_shift);
%     assignin('base', 'V_prof_pre', V_prof_pre);
%     assignin('base', 'cyc_mph', cyc_mph);
%     assignin('base', 'time_final', time_final);
%     
%     Results = sim('ECMS_final','SrcWorkspace', 'base');
%     actual_spd = Results.get('actual_spd');
%     samp_time = evalin('base', 'samp_time');
%     
%     test0_60 = zeros(size(actual_spd))';
%     test0_40 = zeros(size(actual_spd))';
%     test0_85 = zeros(size(actual_spd))';
%     for i= 1:length(actual_spd)
%         if (actual_spd(i)>= 85)
%             test0_85(i) = 1;
%             test0_60(i) = 1;
%             test0_40(i) = 1;
%         elseif (actual_spd(i) >= 60)
%             test0_85(i) = 0;
%             test0_60(i) = 1;
%             test0_40(i) = 1;
%         elseif (actual_spd(i) >= 40)
%             test0_85(i) = 0;
%             test0_60(i) = 0;
%             test0_40(i) = 1;
%         else
%             test0_85(i) = 0;
%             test0_60(i) = 0;
%             test0_40(i) = 0;
%         end
%     end
%     z0_40 = find(test0_40);z0_60 = find(test0_60);z0_85 = find(test0_85);
%     time0_40 = min(z0_40)*samp_time; time0_60 = min(z0_60)*samp_time;
%     time0_85 = min(z0_85)*samp_time; time40_60 = time0_60-time0_40;
%     
% if ~isempty(z0_40)&~isempty(z0_60)&~isempty(z0_85)
%    con(offset+1,1)=time0_60;
%    con(offset+2,1)=time40_60;
%    con(offset+3,1)=time0_85;
% else
%    con(offset+1:offset+3,1)=100;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%    Run Grade Test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     theta = 5;
%     assignin('base', 'theta', theta); 
%     init_soc = 0.65;
%     assignin('base', 'ess_init_soc', init_soc);
%     timetomaintain = 100;   % in sec
%     test_spd = 55;          % in mph
% 
%     total_test_time = 300;  % in sec
%     ramp_up_time = 40;      % in sec
% 
%     cyc_mph = [0:total_test_time; 0:total_test_time]';
%     cyc_mph(:,2) = test_spd;
%     cyc_mph(1,2) = 0;
%     for i=1:(ramp_up_time)
%         cyc_mph(i+1,2) = i*(test_spd/ramp_up_time);
%     end
%     time_final = length(cyc_mph(:,2));
%     V_prof=cyc_mph;
%     time_shift=1;
%     V_prof_pre=[cyc_mph(1+time_shift:end,1); zeros(time_shift,1)];
%     
%     assignin('base', 'V_prof',  V_prof);
%     assignin('base', 'time_shift', time_shift);
%     assignin('base', 'V_prof_pre', V_prof_pre);
%     
%     assignin('base', 'cyc_mph', cyc_mph);
%     assignin('base', 'time_final', time_final);
%     
%     Results = sim('ECMS_final','SrcWorkspace', 'base');   
%     actual_spd = Results.get('actual_spd');
%     
%     counter = 0;
%     pass_test = [0];
%     threshold_spd = 0.95*test_spd;
%     test_start_time = find(actual_spd >= threshold_spd,1);
%     samp_time = evalin('base', 'samp_time');
% if test_start_time ~= 0
%     for i = test_start_time:length(actual_spd)
%         if (actual_spd(i) >= threshold_spd)
%             counter = counter + samp_time;
%         else
%             pass_test(length(pass_test) + 1) = counter;
%             counter = 0;
%         end
%     end
%     pass_test(length(pass_test)+1) = counter;
%     max_time_above = max(pass_test); % in sec
%     
%     if (max_time_above >= timetomaintain)
%         con(offset+4,1)=theta;
%     else
%         con(offset+4,1)=0;
%     end
% else
%     con(offset+4,1)=0;
% end

con_e=0;

return
