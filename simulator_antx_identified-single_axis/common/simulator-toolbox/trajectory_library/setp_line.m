function set_point = setp_line( t, simulation_time, pos_ned_init, attitude_init, des_att, altitude, speed )
%The movement is composed by the following parts:
%1- Quadrotor starts and remains in the initial pos/att
%2- Use of two 5th order polynomials (1 for pos and 1 for att) to move the
%   drone from initial pos/att to [0;0;height] and null attitude
%3- Hovering in [0;0; height] and null attitude
%4- position step from the center of the cage to the desired point

%% Timing and conditions
t1 = 5;                     %Starting time with drone on the floor
t2 = t1 + 10;               %Time to take off and go to [0;0;height] and null att.
t3 = t2 + 5;                %Time for hovering in [0;0;height]
t4 = t3 + 5;                %Time to accelerate
t5 = simulation_time - 10;  %Constant speed time
t6 = simulation_time - 5;   %Time do decelerate and land

%% Position set-point

%Conditions at the beginning of simulation
pos_i = pos_ned_init;       %Initial position
vel_i = zeros(3,1);         %Initial speed
acc_i = zeros(3,1);         %Initial acceleration
att_i = attitude_init;      %Initial attitude
ang_rates_i = zeros(3,1);   %Initial angular speed
ang_acc_i = zeros(3,1);     %Initial angular acceleration

%Conditions of hovering
pos_h = pos_ned_init + [0; 0; altitude];    %Hovering position
vel_h = zeros(3,1);        %Hovering speed
acc_h = zeros(3,1);        %Hovering acceleration
att_h = des_att;           %Hovering attitude
ang_rates_h = zeros(3,1);  %Hovering angular speed
ang_acc_h = zeros(3,1);    %Hovering angular acceleration

%Conditions at the end of acceleration
pos_a = pos_h + speed * (t4 - t3);    %Hovering position
vel_a = speed;                  %Hovering speed
acc_a = zeros(3,1);             %Hovering acceleration
att_a = des_att;                %Hovering attitude
ang_rates_a = zeros(3,1);       %Hovering angular speed
ang_acc_a = zeros(3,1);         %Hovering angular acceleration

%Deceleration conditions
pos_di = pos_a + speed * (t5 - t4);
pos_df = pos_di + [0.1 * pos_di(1,1); 0.1 * pos_di(2,1); 0];
vel_di = speed;
vel_df = zeros(3,1);
acc_di = acc_a;
acc_df = zeros(3,1);
att_di = att_a;
att_df = zeros(3,1);
ang_rates_di = ang_rates_a;
ang_rates_df = zeros(3,1);
ang_acc_di = ang_acc_a;
ang_acc_df = zeros(3,1);



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
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t1, t2, pos_i, pos_h, vel_i, vel_h, acc_i, acc_h);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t1, t2, att_i, att_h, ang_rates_i, ang_rates_h, ang_acc_i, ang_acc_h);
    
elseif t >= t2 && t < t3
    pos_sp = pos_h;
    vel_sp = vel_h;
    acc_sp = acc_h;
    jerk_sp = zeros(3,1);
    attitude_sp = att_h;
    ang_rates_sp = ang_rates_h;
    ang_acc_sp = ang_acc_h;
    ang_jerk_sp = zeros(3,1);
    
elseif t >= t3 && t < t4
    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t3, t4, pos_h, pos_a, vel_h, vel_a, acc_h, acc_a);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t3, t4, att_h, att_a, ang_rates_h, ang_rates_a, ang_acc_h, ang_acc_a);
    
elseif t >= t4 && t < t5
    
    pos_sp = pos_a + speed * (t - t4);
    vel_sp = speed;
    acc_sp = acc_a;
    jerk_sp = zeros(3,1);
    attitude_sp = att_a;
    ang_rates_sp = ang_rates_a;
    ang_acc_sp = ang_acc_a;
    ang_jerk_sp = zeros(3,1);
    
elseif t >= t5 && t < t6
      
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t5, t6, pos_di, pos_df, vel_di, vel_df, acc_di, acc_df);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t5, t6, att_di, att_df, ang_rates_di, ang_rates_df, ang_acc_di, ang_acc_df);
    
else
        
    pos_sp = pos_df;
    vel_sp = zeros(3,1);
    acc_sp = zeros(3,1);
    jerk_sp = zeros(3,1);
    attitude_sp = zeros(3,1);
    ang_rates_sp = zeros(3,1);
    ang_acc_sp = zeros(3,1);
    ang_jerk_sp = zeros(3,1);
    
end

set_point = [pos_sp; vel_sp; acc_sp; jerk_sp; ...
    attitude_sp; ang_rates_sp; ang_acc_sp; ang_jerk_sp];
end