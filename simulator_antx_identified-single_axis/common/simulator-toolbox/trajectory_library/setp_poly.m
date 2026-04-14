function set_point = setp_poly( t, simulation_time, pos_ned_init, pos_ned_fin, attitude_init, height )
%The movement is composed by 9 parts:
%1- Tiltrotor starts and remains in the initial pos/att
%2- Use of two 5th order polynomials (1 for pos and 1 for att) to move the
%   drone from initial pos/att to [0;0;height] and null attitude
%3- Hovering in [0;0; height] and null attitude
%4- Use of a 5th order polynomial in order to move from null attitude to the
%first attitude change: +90° yaw
%5-Use of a 5th order polynomial in order to move from first attitude to the
%seconf attitude change: +45° pitch
%6- hovering with non null attitude
%7- 5th order polys to return with null attitude in [0;0;0]
%8- Landing in [0;0;0] with null attitude
%9- "Hovering" in [0;0;0] with null attitude

%% Timing and conditions
t1 = 5;                     %Starting time with drone on the floor
t2 = t1 + 10;               %Time to take off and go to [0;0;height] and null att.
t3 = t2 + 5;                %Time for hovering in [0;0;height]
t4 = t3 + 5;                %Time to change the attitude staying in [0;0;height]: +90° of yaw
t5 = t4 + 10;               %Time to change attitude: +45° pitch. Then hovering
t6 = simulation_time-25;    %Time to return to null attitude
t7 = simulation_time-15;    %Time to land
t8 = simulation_time-5;     %5 secs hovering on the ground

%% Position set-point

%Conditions at the beginning of simulation
pos_i = pos_ned_init;       %Initial position
vel_i = zeros(3,1);         %Initial speed
acc_i = zeros(3,1);         %Initial acceleration
att_i = attitude_init;      %Initial attitude
ang_rates_i = zeros(3,1);   %Initial angular speed
ang_acc_i = zeros(3,1);     %Initial angular acceleration

%Conditions of hovering at the center of the cage
pos_h = [0; 0; height];    %Hovering position
vel_h = zeros(3,1);        %Hovering speed
acc_h = zeros(3,1);        %Hovering acceleration
att_h = [0; 0; 0];         %Hovering attitude
ang_rates_h = zeros(3,1);  %Hovering angular speed
ang_acc_h = zeros(3,1);    %Hovering angular acceleration

%First attitude variation: +90° yaw
pos_a1 = [0;0;height];
vel_a1 = zeros(3,1);        
acc_a1 = zeros(3,1);        
att_a1 = [0; 0; 90]*pi/180;         
ang_rates_a1 = zeros(3,1);  
ang_acc_a1 = zeros(3,1);    

%Second attitude variation: +45° pitch
pos_a2 = [0;0;height];
vel_a2 = zeros(3,1);        
acc_a2 = zeros(3,1);        
att_a2 = att_a1 + [0; 45; 0]*pi/180;  
ang_rates_a2 = zeros(3,1);  
ang_acc_a2 = zeros(3,1);    

%Landing conditions
pos_l = pos_ned_fin;           
vel_l = zeros(3,1);         
acc_l = zeros(3,1);         
att_l = zeros(3,1);          
ang_rates_l = zeros(3,1);   
ang_acc_l = zeros(3,1);     


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
    %first attitude change
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t3, t4, pos_h, pos_a1, vel_h, vel_a1, acc_h, acc_a1);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t3, t4, att_h, att_a1, ang_rates_h, ang_rates_a1, ang_acc_h, ang_acc_a1);

elseif t >= t4 && t < t5
    %second attitude change
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t4, t5, pos_a1, pos_a2, vel_a1, vel_a2, acc_a1, acc_a2);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t4, t5, att_a1, att_a2, ang_rates_a1, ang_rates_a2, ang_acc_a1, ang_acc_a2);

elseif t >= t5 && t < t6
    pos_sp = pos_a2;
    vel_sp = vel_a2;
    acc_sp = acc_a2;
    jerk_sp = zeros(3,1);
    attitude_sp = att_a2; 
    ang_rates_sp = ang_rates_a2; 
    ang_acc_sp = ang_acc_a2; 
    ang_jerk_sp = zeros(3,1);

elseif t >= t6 && t < t7
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t6, t7, pos_a2, pos_h, vel_a2, vel_h, acc_a2, acc_h);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t6, t7, att_a2, att_h, ang_rates_a2, ang_rates_h, ang_acc_a2, ang_acc_h);

elseif t >= t7 && t < t8    
    [pos_sp, vel_sp, acc_sp, jerk_sp] = fifth_ord_poly(t, t7, t8, pos_h, pos_l, vel_h, vel_l, acc_h, acc_l);
    [attitude_sp, ang_rates_sp, ang_acc_sp, ang_jerk_sp] = fifth_ord_poly(t, t7, t8, att_h, att_l, ang_rates_h, ang_rates_l, ang_acc_h, ang_acc_l);

else    
    pos_sp = pos_l;
    vel_sp = vel_l;
    acc_sp = acc_l;
    jerk_sp = zeros(3,1);
    attitude_sp = att_l; 
    ang_rates_sp = ang_rates_l; 
    ang_acc_sp = ang_acc_l; 
    ang_jerk_sp = zeros(3,1);

end

set_point = [pos_sp; vel_sp; acc_sp; jerk_sp; ...
    attitude_sp; ang_rates_sp; ang_acc_sp; ang_jerk_sp];
end