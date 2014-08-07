% New for this Release

% Made better plots

% Found a major bug --- fixed it with:

                    % Update x1
                    eff_m = interp2(m_map_trq, m_map_spd, m_eff_map, Tm_c, abs(Wm_c))';      
                    eff_m(isnan(eff_m)) = 0.2;
                    Pm_c = (Wm_c*Tm_c).*(eff_m.^-sign(Tm_c));                             % If the torque is negative, we are charging the battery
                    Pbat = repmat(Pm_c,[1, x1_length]);                                      
                    Pbat = permute(Pbat, [2 1]);
                    
                    % Discharge
                    Pbatt_max = repmat(interp1(ess_soc, ess_max_pwr_dis, x1_grid),[1,u1_length]);
                    infeasible_Pbatt(:,x2,x3,:,u2,u3) = (Pbatt_max < Pbat);
                    Pbat(Pbatt_max < Pbat) = Pbatt_max(Pbatt_max < Pbat);
                    rint_discharge = repmat(interp1(ess_soc,ess_r_dis,x1_grid),[1,u1_length]);
                    
                    % Charge
                    Pbatt_min = -repmat(interp1(ess_soc, ess_max_pwr_chg, x1_grid),[1,u1_length]);
                    infeasible_Pbatt(:,x2,x3,:,u2,u3) = (Pbat < Pbatt_min);    % Can brake to make up the difference
                    Pbat(Pbat < Pbatt_min) = Pbatt_min(Pbat < Pbatt_min);
                    rint_charge = repmat(interp1(ess_soc,ess_r_chg,x1_grid),[1,u1_length]);
                  
                    % Charge & Discharge Resistances
                    rint_c = rint_charge;
                    rint_c(Pbat > 0) = rint_discharge(Pbat > 0);
                    
                    
          %   --- It was always going into the else statement because now
          %   all of Pbatt was positive!! - was not saturating it with Pbat
          %   max and then we were getting negative SOC's!!
          
          
          % THis also fixed the SOC, so that there is no more imaginary