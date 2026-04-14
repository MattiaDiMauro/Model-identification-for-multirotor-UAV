function setpoint = setp_spiral( t, simulation_time, pos_ned_init, pos_ned_fin, attitude_init, altitude, omega, roll_pitch, yaw_following, scale )

% [ pos_d, vel_d, acc_d, jerk_d ] =  fifth_ord_poly ( t, ti, tf, xi, xf, vi, vf, ai, af )

%% Timing and conditions

t1 = 5;                     %Starting time with drone on the floor
t2 = t1 + 5;                %Time to take off and go to [0;0;height] and null att.
t3 = t2 + 5;                %Time for hovering in [0;0;height]
t4 = t3 + 5 * (1.4 - omega);%Time for moving from hovering point to the beginning of the spiral traj.
t5 = simulation_time-15;    %Spiral trajectory
t6 = simulation_time-8;    %Hovering and then start to land
t7 = simulation_time-4;     %End time of landing procedure
%At the end 5 seconds of hovering in the ground

b=-0.05;                    %Convergence speed
h = 1.5;                    %Altitude of the christmas tree

%Conditions at the beginning of the simulation
pos_i = pos_ned_init;       %Initial position
vel_i = zeros(3,1);         %Initial speed
acc_i = zeros(3,1);         %Initial acceleration
att_i = attitude_init;      %Initial attitude
ang_rates_i = zeros(3,1);   %Initial angular speed
ang_acc_i = zeros(3,1);     %Initial angular acceleration

%Conditions of the hovering at the center of the cage
pos_h = [0; 0; altitude];   %Hovering position
vel_h = zeros(3,1);         %Hovering speed
acc_h = zeros(3,1);         %Hovering acceleration
att_h = [0; 0; 0];          %Hovering attitude
ang_rates_h = zeros(3,1);   %Hovering angular speed
ang_acc_h = zeros(3,1);     %Hovering angular acceleration

%Conditions at the beginning the spiral trajectory
pos_spiral_i = scale * [1; 0; altitude / scale];                    %Initial position of spiral trajectory
vel_spiral_i = scale * [b; omega; 0];                               %Initial velocity of spiral trajectory
acc_spiral_i = scale * [(b + omega)*(b - omega); 2*b*omega; 0];     %Initial acceleraton of spiral trajectory
att_spiral_i = [-roll_pitch; roll_pitch; -pi/4];                    %Initial attitude of spiral trajectory
if yaw_following == 1
    ang_rates_spiral_i = [0; 0; omega];        %Initial angular speed of spiral trajectory
    ang_acc_spiral_i = [0; 0; 0];                            %Initial angular acceleration of spiral trajectory
else
    ang_rates_spiral_i = [0; 0; 0];             %Initial angular speed of spiral trajectory
    ang_acc_spiral_i = [0; 0; 0];                            %Initial angular acceleration of spiral trajectory
end

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
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t3, t4, pos_h, pos_spiral_i, vel_h, vel_spiral_i, acc_h, acc_spiral_i);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t3, t4, att_h, att_spiral_i, ang_rates_h, ang_rates_spiral_i, ang_acc_h, ang_acc_spiral_i);
    
elseif (t >= t4 && t < t5)
    pos_sp = scale * [ exp( b * (t - t4) ) * cos( omega * (t - t4) );
        exp( b * (t - t4) ) * sin( omega * (t - t4) );
        altitude / scale             ];
    vel_sp = scale * [ exp( b * (t - t4) ) * (b * cos( omega * (t - t4) ) - omega * sin( omega * (t - t4) ));
        exp( b * (t - t4) ) * (b * sin( omega * (t - t4) ) + omega * cos( omega * (t - t4) ));
        0                                                            ];
    acc_sp = scale * [ exp( b * (t - t4)) * (( b - omega ) * (b + omega) * cos ( omega * (t - t4) ) - 2 * b * omega * sin( omega * (t - t4) ));
        exp( b * (t - t4)) * (( b - omega ) * (b + omega) * sin ( omega * (t - t4) ) + 2 * b * omega * cos( omega * (t - t4) ));
        0                                                                                              ];
    jerk_sp = scale * [ exp( b * (t - t4) ) * ( omega * ( omega^2 - 3 * b^2 ) * sin( omega * (t - t4) ) + b * (b^2 - 3 * omega^2) * cos(omega * (t - t4)));
        exp( b * (t - t4) ) * ( b * ( b^2 - 3 * omega^2 ) * sin( omega * (t - t4) ) - omega * (omega^2 - 3 * b^2) * cos(omega * (t - t4)));
        0                                                                                                          ];
    
    [alt, alt_dot, alt_ddot, alt_dddot] = fifth_ord_poly(t, t4, t5, [0;0;0], [0;0;-h], [0;0;0], [0;0;0], [0;0;0], [0;0;0]);
    
    pos_sp(3) = pos_sp(3) + alt(3);
    vel_sp(3) = vel_sp(3) + alt_dot(3);
    acc_sp(3) = acc_sp(3) + alt_ddot(3);
    jerk_sp(3) = jerk_sp(3) + alt_dddot(3);
    
    if yaw_following == 1
        attitude_sp =   [-roll_pitch; roll_pitch; omega * (t - t4) - pi/4];
        ang_rates_sp =  [0; 0; omega];
        ang_acc_sp =    [0; 0; 0];
        ang_jerk_sp =   [0; 0; 0];
    else
        attitude_sp =   [-roll_pitch; roll_pitch; - pi/4];
        ang_rates_sp =  [0; 0; 0];
        ang_acc_sp =    [0; 0; 0];
        ang_jerk_sp =   [0; 0; 0];
    end
    
elseif (t >= t5 && t < t6)
    
    pos_spiral_f = scale * [ exp( b * (t5 - t4) ) * cos( omega * (t5 - t4) );
        exp( b * (t5 - t4) ) * sin( omega * (t5 - t4) );
        altitude / scale - h / scale            ];
    vel_spiral_f = scale * [ exp( b * (t5 - t4) ) * (b * cos( omega * (t5 - t4) ) - omega * sin( omega * (t5 - t4) ));
        exp( b * (t5 - t4) ) * (b * sin( omega * (t5 - t4) ) + omega * cos( omega * (t5 - t4) ));
        0                                                                  ];
    acc_spiral_f = scale * [ exp( b * (t5 - t4)) * (( b - omega ) * (b + omega) * cos ( omega * (t5 - t4) ) - 2 * b * omega * sin( omega * (t5 - t4) ));
        exp( b * (t5 - t4)) * (( b - omega ) * (b + omega) * sin ( omega * (t5 - t4) ) + 2 * b * omega * cos( omega * (t5 - t4) ));
        0                                                                                                       ];
    jerk_spiral_f = scale * [ exp( b * (t5 - t4) ) * ( omega * ( omega^2 - 3 * b^2 ) * sin( omega * (t5 - t4) ) + b * (b^2 - 3 * omega^2) * cos(omega * (t5 - t4)));
        exp( b * (t5 - t4) ) * ( b * ( b^2 - 3 * omega^2 ) * sin( omega * (t5 - t4) ) - omega * (omega^2 - 3 * b^2) * cos(omega * (t5 - t4)));
        0                                    ];
    
    if yaw_following == 1
        att_spiral_f =         [-roll_pitch; roll_pitch; omega*(t5-t4) - floor(omega*(t5-t4)/4/pi)*4*pi - pi/4];
        ang_rates_spiral_f =   [0; 0; omega ];
        ang_acc_spiral_f =     [0; 0; 0];
        
        if mod(floor(omega*(t5-t4)/2/pi),2)
            att_h(3) = att_h(3)+2*pi;
        end
        
    else
        att_spiral_f =         [-roll_pitch; roll_pitch; - pi/4];
        ang_rates_spiral_f =   [0; 0; 0];
        ang_acc_spiral_f =     [0; 0; 0];
    end
    
    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t5, t6, pos_spiral_f, pos_h, vel_spiral_f, vel_h, acc_spiral_f, acc_h);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t5, t6, att_spiral_f, att_h, ang_rates_spiral_f, ang_rates_h, ang_acc_spiral_f, ang_acc_h);
    
elseif (t >= t6 && t < t7)
    if yaw_following == 1
        if mod(floor((omega*(t5-t4))/pi),4) == 1 || mod(floor((omega*(t5-t4))/pi),4) == 2 
            att_h(3) = att_h(3)+2*pi;
            att_l(3) = att_l(3) + 2*pi;
        end
    end
    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t6, t7, pos_h, pos_l, vel_h, vel_l, acc_h, acc_l);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t6, t7, att_h, att_l, ang_rates_h, ang_rates_l, ang_acc_h, ang_acc_l);
    
else
    if yaw_following == 1
        if mod(floor((omega*(t5-t4))/pi),4) == 1 || mod(floor((omega*(t5-t4))/pi),4) == 2
            att_l(3) = att_l(3) + 2*pi;
        end
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

