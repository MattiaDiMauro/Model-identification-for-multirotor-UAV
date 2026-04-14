function [ omega, omega_dot ] = fromNEDtoBODYangvel( attitude, attitude_dot, attitude_2dot )

phi = attitude(1);
theta = attitude(2);

phi_dot = attitude_dot(1);
theta_dot = attitude_dot(2);

G = [1 0 -sin(theta); 
     0 cos(phi) sin(phi)*cos(theta);
     0 -sin(phi) cos(phi)*cos(theta)];

G_dot = [0 0 -theta_dot*cos(theta); 
         0 -phi_dot*sin(phi) phi_dot*cos(phi)*cos(theta)-theta_dot*sin(phi)*sin(theta); 
         0 -phi_dot*cos(phi) -phi_dot*sin(phi)*cos(theta)-theta_dot*sin(theta)*cos(phi)];

omega = G * attitude_dot;     
omega_dot = G * attitude_2dot + G_dot * attitude_dot;

end

