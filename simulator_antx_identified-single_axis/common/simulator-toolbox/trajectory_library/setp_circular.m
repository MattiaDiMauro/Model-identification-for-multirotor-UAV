function setpoint = setp_circular( t, simulation_time, pos_ned_init, pos_ned_fin, attitude_init, altitude, omega, roll, yaw_following, scale )

%The movement is composed by 7 parts:
%1- Tiltrotor starts and remains in the initial pos/att
%2- Use of two 5th order polynomials (1 for pos and 1 for att) to move the
%   drone from initial pos/att to [0;0;height] and null attitude
%3- Hovering in [0;0; height] and null attitude
%4- Use of two 5th order polynomial in order to move from [0;0;height] to the
%beginning of the circular trajectory.
%5- circular trajectory with roll movement
%6- landing in [0;0;0] with null attitude
%7- "Hovering" in [0;0;0] with null attitude

%IF yaw_following == 2, the yaw angle would follow the tangent direction of the circular trajectory
% and the roll angle is kept constant at the 'roll' input value

%% Timing and conditions

t1 = 5;                     %Starting time with drone on the floor
t2 = t1 + 10;               %Time to take off and go to [0;0;height] and null att.
t3 = t2 + 5;                %Time for hovering in [0;0;height]
t4 = t3 + 5 * (1.4 - omega);%Time for moving from hovering point to the beginning of the circular traj.
t5 = simulation_time-15;    %End time of 8-shape trajectory
t6 = t5 + 5 * (1.4 - omega);%End time of landing procedure
%At the end 5 seconds of hovering in the ground

%Conditions at the beginning of simulation
pos_i = pos_ned_init;       %Initial position
vel_i = zeros(3,1);         %Initial speed
acc_i = zeros(3,1);         %Initial acceleration
att_i = attitude_init;      %Initial attitude
ang_rates_i = zeros(3,1);   %Initial angular speed
ang_acc_i = zeros(3,1);     %Initial angular acceleration

%Conditions of the hovering at the center of the cage
pos_h = [0; 0; altitude];    %Hovering position
vel_h = zeros(3,1);        %Hovering speed
acc_h = zeros(3,1);        %Hovering acceleration
att_h = [0; 0; 0];         %Hovering attitude
ang_rates_h = zeros(3,1);  %Hovering angular speed
ang_acc_h = zeros(3,1);    %Hovering angular acceleration

%Conditions at the beginning of the circular trajectory
pos_circ_i = scale * [1; 0; altitude / scale];                    %Initial position of circ trajectory
vel_circ_i = scale * [0; omega; 0];                     %Initial speed of circ trajectory
acc_circ_i = scale * [-omega^2 ; 0; 0];                 %Initial acc of circ trajectory
if yaw_following == 1
    att_circ_i = [0; 0; pi/2];                         %Initial att of circ trajectory
    ang_rates_circ_i = [omega * roll; 0; omega];       %Initial angular speed of 8-shape trajectory
    
elseif yaw_following == 2
    att_circ_i = [roll; 0; pi/2];                      %Initial att of circ trajectory
    ang_rates_circ_i = [0; 0; omega];                  %Initial angular speed of 8-shape trajectory
else
    att_circ_i = [0; 0; 0];
    ang_rates_circ_i = [omega * roll; 0; 0];        %Initial angular speed of 8-shape trajectory
end
ang_acc_circ_i = [0;0;0];                       %Initial angular acceleration of 8-shape trajectory

%Landing conditions
pos_l = pos_ned_fin;                    %Landing position
vel_l = zeros(3,1);                     %Landing velocity
acc_l = zeros(3,1);                     %Landing acceleration
att_l = [0;0;0]*pi/180;                 %Landing attitude
ang_rates_l = zeros(3,1);               %Landing angular speed
ang_acc_l = zeros(3,1);                 %Landing angular acceleration

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
    
elseif (t >= t1 && t < t2)
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t1, t2, pos_i, pos_h, vel_i, vel_h, acc_i, acc_h);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t1, t2, att_i, att_h, ang_rates_i, ang_rates_h, ang_acc_i, ang_acc_h);
    
elseif (t >= t2 && t < t3)
    pos_sp = pos_h;
    vel_sp = vel_h;
    acc_sp = acc_h;
    jerk_sp = zeros(3,1);
    attitude_sp = att_h;
    ang_rates_sp = ang_rates_h;
    ang_acc_sp = ang_acc_h;
    ang_jerk_sp = zeros(3,1);
    
elseif (t >= t3 && t < t4)
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t3, t4, pos_h, pos_circ_i, vel_h, vel_circ_i, acc_h, acc_circ_i);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t3, t4, att_h, att_circ_i, ang_rates_h, ang_rates_circ_i, ang_acc_h, ang_acc_circ_i);
    
elseif (t >= t4 && t < t5)
    pos_sp = scale * [          cos(omega * (t - t4)) ;
        sin(omega * (t - t4)) ;
        altitude / scale       ];
    vel_sp = scale * [   -omega * sin(omega * (t - t4)) ;
        omega * cos(omega * (t - t4)) ;
        0  ];
    acc_sp = scale * [        -omega^2 * cos(omega * (t - t4)) ;
        -omega^2 * sin(omega * (t - t4)) ;
        0  ];
    jerk_sp = scale * [         omega^3 * sin(omega * (t - t4)) ;
        -omega^3 * cos(omega * (t - t4)) ;
        0  ];
    if yaw_following == 1
        attitude_sp = [roll*sin(omega*(t-t4)); 0; pi/2 + omega * (t - t4)];
        ang_rates_sp = [roll*omega*cos(omega*(t-t4)); 0; 1 / ( vel_sp(2)^2 + vel_sp(1)^2 ) * ( acc_sp(2) * vel_sp(1) - vel_sp(2) * acc_sp(1))];
        ang_acc_sp = [-roll*omega^2*sin(omega*(t-t4)); 0; (jerk_sp(2)*vel_sp(1)-vel_sp(2)*jerk_sp(1)) / (vel_sp(1)^2+vel_sp(2)^2) - (acc_sp(2)*vel_sp(1)-vel_sp(2)*acc_sp(1)) * (2*vel_sp(1)*acc_sp(1)+2*vel_sp(2)*acc_sp(2)) / (vel_sp(1)^2+vel_sp(2)^2)^2];
        ang_jerk_sp = zeros(3,1);
        
    elseif yaw_following == 2
        attitude_sp = [roll; 0; pi/2 + omega * (t - t4)];
        ang_rates_sp = [0; 0; 1 / ( vel_sp(2)^2 + vel_sp(1)^2 ) * ( acc_sp(2) * vel_sp(1) - vel_sp(2) * acc_sp(1))];
        ang_acc_sp = [0; 0; (jerk_sp(2)*vel_sp(1)-vel_sp(2)*jerk_sp(1)) / (vel_sp(1)^2+vel_sp(2)^2) - (acc_sp(2)*vel_sp(1)-vel_sp(2)*acc_sp(1)) * (2*vel_sp(1)*acc_sp(1)+2*vel_sp(2)*acc_sp(2)) / (vel_sp(1)^2+vel_sp(2)^2)^2];
        ang_jerk_sp = zeros(3,1);
    else
        attitude_sp = [roll*sin(omega*(t-t4)); 0; 0];
        ang_rates_sp = [roll*omega*cos(omega*(t-t4)); 0; 0];
        ang_acc_sp = [-roll*omega^2*sin(omega*(t-t4)); 0; 0];
        ang_jerk_sp = zeros(3,1);
    end
    
elseif (t >= t5 && t < t6)
    pos_circ_f = scale * [cos(omega * (t5 - t4)); sin( omega * (t5 - t4)); altitude / scale];
    vel_circ_f = scale * [-omega * sin(omega * (t5 - t4)); omega * cos( omega * (t5 - t4)); 0];
    acc_circ_f = scale * [-omega^2 * cos(omega * (t5 - t4));  - omega^2 * sin( omega * (t5 - t4)); 0];
    jerk_circ_f = scale * [omega^3 * sin(omega * (t5 - t4));  - omega^3 * cos( omega * (t5 - t4)); 0];
    if yaw_following == 1
        att_circ_f = [roll*sin(omega*(t5-t4)); 0; atan2(vel_circ_f(2),vel_circ_f(1))];
        ang_rates_circ_f = [roll*omega*cos(omega*(t5-t4)); 0; 1 / ( vel_circ_f(2)^2 + vel_circ_f(1)^2 ) * ( acc_circ_f(2) * vel_circ_f(1) - vel_circ_f(2) * acc_circ_f(1))];
        ang_acc_circ_f = [-roll*omega^2*sin(omega*(t5-t4)); 0; (jerk_circ_f(2)*vel_circ_f(1)-vel_circ_f(2)*jerk_circ_f(1)) / (vel_circ_f(1)^2+vel_circ_f(2)^2) - (acc_circ_f(2)*vel_circ_f(1)-vel_circ_f(2)*acc_circ_f(1)) * (2*vel_circ_f(1)*acc_circ_f(1)+2*vel_circ_f(2)*acc_circ_f(2)) / (vel_circ_f(1)^2+vel_circ_f(2)^2)^2];
    elseif yaw_following == 2
        att_circ_f = [roll; 0; atan2(vel_circ_f(2),vel_circ_f(1))];
        ang_rates_circ_f = [0; 0; 1 / ( vel_circ_f(2)^2 + vel_circ_f(1)^2 ) * ( acc_circ_f(2) * vel_circ_f(1) - vel_circ_f(2) * acc_circ_f(1))];
        ang_acc_circ_f = [0; 0; (jerk_circ_f(2)*vel_circ_f(1)-vel_circ_f(2)*jerk_circ_f(1)) / (vel_circ_f(1)^2+vel_circ_f(2)^2) - (acc_circ_f(2)*vel_circ_f(1)-vel_circ_f(2)*acc_circ_f(1)) * (2*vel_circ_f(1)*acc_circ_f(1)+2*vel_circ_f(2)*acc_circ_f(2)) / (vel_circ_f(1)^2+vel_circ_f(2)^2)^2];
    else
        att_circ_f = [ roll*sin(omega*(t5-t4)); 0; 0 ];
        ang_rates_circ_f = [ roll*omega*cos(omega*(t5-t4)); 0; 0 ];
        ang_acc_circ_f = [-roll*omega^2*sin(omega*(t5-t4)); 0; 0];
    end
    
    if (yaw_following == 1 || yaw_following == 2) && (mod(floor((pi/2+omega*(t5-t4))/pi),4) == 1 || mod(floor((pi/2+omega*(t5-t4))/pi),4) == 2)
        att_circ_f = [att_circ_f(1); att_circ_f(2); att_circ_f(3)+2*pi];
        att_l = [att_l(1); att_l(2); att_l(3)+2*pi];
    end
    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t5, t6, pos_circ_f, pos_l, vel_circ_f, vel_l, acc_circ_f, acc_l);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t5, t6, att_circ_f, att_l, ang_rates_circ_f, ang_rates_l, ang_acc_circ_f, ang_acc_l);
    
else
    
    if (yaw_following == 1 || yaw_following == 2) && (mod(floor((pi/2+omega*(t5-t4))/pi),4) == 1 || mod(floor((pi/2+omega*(t5-t4))/pi),4) == 2)
        att_l = [att_l(1); att_l(2); att_l(3)+2*pi];
    end
    
    pos_sp = pos_l;
    vel_sp = vel_l;
    acc_sp = acc_l;
    jerk_sp = zeros(3,1);
    attitude_sp =  att_l;
    ang_rates_sp =  ang_rates_l;
    ang_acc_sp = ang_acc_l;
    ang_jerk_sp = zeros(3,1);
    
end

%% Setpoint output

setpoint = [pos_sp; vel_sp; acc_sp; jerk_sp; ...
    attitude_sp; ang_rates_sp; ang_acc_sp; ang_jerk_sp];

end

