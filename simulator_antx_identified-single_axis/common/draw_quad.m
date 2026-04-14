function draw_quad(input, arm_length, propeller_radius)
% ----------------------------------------------------
% QUADROTOR VERSION
% Draws the quadrotor in 3 dimensions during the simulation.
% This script is invoked once for each animation frame.
% ----------------------------------------------------
%#codegen

global old_position;
global index_view;
global quadrotor;

THROTTLE_MAX = 0.3;

% NED position and attitude
n = input(1);
e = input(2);
d = input(3);
quatNed = input(4:7);

% Convert position and attitude from NED to ENU (a.k.a. XYZ)
x = e;
y = n;
z = -d;
phi = pi;
qx = [sin(phi/2) ;
          0      ;
          0      ;
      cos(phi/2)];

psi = -pi/2;
qz = [    0      ;
          0      ;
      sin(psi/2) ;
      cos(psi/2)];

q_enu_to_ned = quatProd( qz, qx );

quat = quatProd( quatNed, quatConj( q_enu_to_ned ) );
norm = sqrt(quat(1)^2 + quat(2)^2 + quat(3)^2 + quat(4)^2);
quatEnu = quat / norm;

% Convert quaternion to attitude matrix
rotation = quatToAtt(quatEnu)';

% Actuators
th_1 = input(8);
th_2 = input(9);
th_3 = input(10);
th_4 = input(11);

if index_view == 0% Code executed only the first time that the script is called
    %% Draw fixed frame reference
    
    %N - Axes
    line(   [0, 0  ], ...
            [0, 0.5], ...
            [0, 0  ], ...
        'linewidth', 2, 'color', 'red');
    text(0, 0.6, 0, 'N', 'fontsize', 13);
    
    %E - Axes
    line(   [0, 0.5], ...
        	[0, 0  ], ...
            [0, 0  ], ...
        'linewidth', 2, 'color', 'red');
    text(0.6, 0, 0, 'E', 'fontsize', 13);
    
    %D - Axes
    line(   [0, 0  ], ...
            [0, 0  ], ...
            [0, -0.5], ...
        'linewidth', 2, 'color', 'red');
    text(0, 0, -0.6, 'D', 'fontsize', 13);
    
    %Ground square
    line(   [-6, 6], ...
            [ 3, 3], ...
            [ 0, 0], ...
        'linewidth',2,'color','black');
    line(   [-6, 6], ...
            [-3,-3], ...
            [   0,   0], ...
        'linewidth',2,'color','black');
    line(   [ 6, 6], ...
            [-3, 3], ...
            [ 0, 0], ...
        'linewidth',2,'color','black');
    line(   [-6,-6], ...
            [-3, 3], ...
            [ 0, 0], ...
        'linewidth',2,'color','black');
    
    % Set the camera position and target
    camtarget_x = 0;
    camtarget_y = 0;
    camtarget_z = 2;
    campos_x = -20;
    campos_y =  -62;
    campos_z = 26;
    
    camtarget([camtarget_x,camtarget_y,camtarget_z]);
    campos([campos_x,campos_y,campos_z]);
    
else % This part is not executed the first time the script is called.
    %% Delete the quadrotor drawing in the old position
    drawnow;
    axis([-6 6 -3 3 -0.6 4]);
    
    delete(quadrotor.arm1);
    delete(quadrotor.arm2);
    delete(quadrotor.arm3);
    delete(quadrotor.arm4);
    delete(quadrotor.prop1);
    delete(quadrotor.prop2);
    delete(quadrotor.prop3);
    delete(quadrotor.prop4);
    delete(quadrotor.th1);
    delete(quadrotor.th2);
    delete(quadrotor.th3);
    delete(quadrotor.th4);
    
    line(   [old_position(1), x], ...
            [old_position(2), y], ...
            [old_position(3), z], ...
        'linewidth',1,'color','magenta');
end

%% Draw the quadrotor frame
arm_1 = [   0, +arm_length * sqrt(2) / 2 ;
            0, +arm_length * sqrt(2) / 2 ;
            0,                         0];
        
arm_2 = [   0, -arm_length * sqrt(2) / 2 ;
            0, -arm_length * sqrt(2) / 2 ;
            0,                         0];
                
arm_3 = [   0, +arm_length * sqrt(2) / 2 ;
            0, -arm_length * sqrt(2) / 2 ;
            0,                         0];

arm_4 = [   0, -arm_length * sqrt(2) / 2 ;
            0, +arm_length * sqrt(2) / 2 ;
            0,                         0];
        
quadrotor.arm1 = line( x + rotation(1,:) * arm_1, ...
                       y + rotation(2,:) * arm_1, ...
                       z + rotation(3,:) * arm_1, ...
                       'linewidth',2,'color','black');

quadrotor.arm2 = line( x + rotation(1,:) * arm_2, ...
                       y + rotation(2,:) * arm_2, ...
                       z + rotation(3,:) * arm_2, ...
                       'linewidth',2,'color','black');

quadrotor.arm3 = line( x + rotation(1,:) * arm_3, ...
                       y + rotation(2,:) * arm_3, ...
                       z + rotation(3,:) * arm_3, ...
                       'linewidth',2,'color','black');

quadrotor.arm4 = line( x + rotation(1,:) * arm_4, ...
                       y + rotation(2,:) * arm_4, ...
                       z + rotation(3,:) * arm_4, ...
                       'linewidth',2,'color','black');

%% Propellers
theta = 0: pi/6 :2*pi;

circ_x = propeller_radius * cos(theta);
circ_y = propeller_radius * sin(theta);
circ_z = zeros(size(theta));

points = [  circ_x ;
            circ_y ;
            circ_z];
                        
quadrotor.prop1 = line( x + rotation(1,:) * arm_1(:,2) + rotation(1,:) * points, ...
                        y + rotation(2,:) * arm_1(:,2) + rotation(2,:) * points, ...
                        z + rotation(3,:) * arm_1(:,2) + rotation(3,:) * points, ...
                        'linewidth',2,'color','red');

quadrotor.prop2 = line( x + rotation(1,:) * arm_2(:,2) + rotation(1,:) * points, ...
                        y + rotation(2,:) * arm_2(:,2) + rotation(2,:) * points, ...
                        z + rotation(3,:) * arm_2(:,2) + rotation(3,:) * points, ...
                        'linewidth',2,'color','blue');

quadrotor.prop3 = line( x + rotation(1,:) * arm_3(:,2) + rotation(1,:) * points, ...
                        y + rotation(2,:) * arm_3(:,2) + rotation(2,:) * points, ...
                        z + rotation(3,:) * arm_3(:,2) + rotation(3,:) * points, ...
                        'linewidth',2,'color','red');

quadrotor.prop4 = line( x + rotation(1,:) * arm_4(:,2) + rotation(1,:) * points, ...
                        y + rotation(2,:) * arm_4(:,2) + rotation(2,:) * points, ...
                        z + rotation(3,:) * arm_4(:,2) + rotation(3,:) * points, ...
                        'linewidth',2,'color','blue');

%% Thrust
thrust_1 = [0,                        0 ;
            0,                        0 ;
            0, -THROTTLE_MAX * th_1/100];

thrust_2 = [0,                        0 ;
            0,                        0 ;
            0, -THROTTLE_MAX * th_2/100];
        
thrust_3 = [0,                        0 ;
            0,                        0 ;
            0, -THROTTLE_MAX * th_3/100];
        
thrust_4 = [0,                        0 ;
            0,                        0 ;
            0, -THROTTLE_MAX * th_4/100];
        
quadrotor.th1 = line( x + rotation(1,:) * arm_1(:,2) + rotation(1,:) * thrust_1, ...
                      y + rotation(2,:) * arm_1(:,2) + rotation(2,:) * thrust_1, ...
                      z + rotation(3,:) * arm_1(:,2) + rotation(3,:) * thrust_1, ...
                      'linewidth',2,'color','red');

quadrotor.th2 = line( x + rotation(1,:) * arm_2(:,2) + rotation(1,:) * thrust_2, ...
                      y + rotation(2,:) * arm_2(:,2) + rotation(2,:) * thrust_2, ...
                      z + rotation(3,:) * arm_2(:,2) + rotation(3,:) * thrust_2, ...
                      'linewidth',2,'color','blue');

quadrotor.th3 = line( x + rotation(1,:) * arm_3(:,2) + rotation(1,:) * thrust_3, ...
                      y + rotation(2,:) * arm_3(:,2) + rotation(2,:) * thrust_3, ...
                      z + rotation(3,:) * arm_3(:,2) + rotation(3,:) * thrust_3, ...
                      'linewidth',2,'color','red');

quadrotor.th4 = line( x + rotation(1,:) * arm_4(:,2) + rotation(1,:) * thrust_4, ...
                      y + rotation(2,:) * arm_4(:,2) + rotation(2,:) * thrust_4, ...
                      z + rotation(3,:) * arm_4(:,2) + rotation(3,:) * thrust_4, ...
                      'linewidth',2,'color','blue');

%% Save the current position for the path plot
old_position = [x, y, z];

% Count the iterations
index_view = index_view + 1;