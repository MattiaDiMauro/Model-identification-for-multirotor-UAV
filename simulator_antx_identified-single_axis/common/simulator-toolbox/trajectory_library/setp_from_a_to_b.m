function set_point = setp_from_a_to_b(t, init_pos, fin_pos, init_att, fin_att, simulation_time)
%setpoint to move from point A to point B, given initial and final
%position and orientation 

%% Timing and conditions
                                                    
t1 = 5;
t2 = simulation_time - 5 ;
%% Position set-point

%Conditions at the beginning of simulation

pos_i = init_pos;           %Initial position
vel_i = zeros(3,1);         %Initial speed
acc_i = zeros(3,1);         %Initial acceleration
att_i = init_att;           %Initial attitude
ang_rates_i = zeros(3,1);   %Initial angular speed
ang_acc_i = zeros(3,1);     %Initial angular acceleration

%Conditions at the end of the simulation

pos_f = fin_pos;           %Final position
vel_f = zeros(3,1);        %Final speed
acc_f = zeros(3,1);        %Final acceleration
att_f = fin_att;           %Final attitude
ang_rates_f = zeros(3,1);  %Final angular speed
ang_acc_f = zeros(3,1);    %Final angular acceleration


%% Trajectory generation

if t < t1
    pos_sp = pos_i;
    vel_sp = vel_i;
    acc_sp = acc_i;
    jerk_sp = zeros(3,1);
    attitude_sp = att_i; 
    ang_rates_sp = ang_rates_i; 
    ang_acc_sp = ang_acc_i; 
    ang_jerk_sp = zeros(3,1);
    
    
elseif t >= t1 && t < t2
    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t1, t2, pos_i, pos_f, vel_i, vel_f, acc_i, acc_f);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t1, t2, att_i, att_f, ang_rates_i, ang_rates_f, ang_acc_i, ang_acc_f);

else
    
    pos_sp = pos_f;
    vel_sp = vel_f;
    acc_sp = acc_f;
    jerk_sp = zeros(3,1);
    attitude_sp = att_f; 
    ang_rates_sp = ang_rates_f; 
    ang_acc_sp = ang_acc_f; 
    ang_jerk_sp = zeros(3,1);

end

set_point = [pos_sp; vel_sp; acc_sp; jerk_sp; ...
    attitude_sp; ang_rates_sp; ang_acc_sp; ang_jerk_sp];
end