%% Extract output
sim_out.time = state_output.ang_rates_body.Time;
sim_out.pos_sp = setpoint.Data(:,1:3);
sim_out.vel_sp = setpoint_velocity.Data;
sim_out.att_sp = setpoint_attitude.Data;       %quaternion
sim_out.ang_vel_sp = setpoint_rates.Data;

%clear setpoint setpoint_attitude setpoint_rates setpoint_velocity

sim_out.body_force = body_force.Data;
%clear body_force;

sim_out.inertial_force = inertial_force.Data;
%clear inertial_force;

sim_out.motors = quad_input.Data(:,1:4);
%clear quad_input;

sim_out.pos_ned = state_output.pos_ned.Data;
sim_out.vel_ned = state_output.vel_ned.Data;
sim_out.ang_rates_body = state_output.ang_rates_body.Data;
sim_out.attitude = state_output.attitude.Data;      %quaternion
%clear state_output

sim_out.acc_ned = state_derivative_output.acc_ned.Data;
sim_out.ang_acc_ned = state_derivative_output.ang_acc_ned.Data;
%clear state_derivative_output

% Convert quaternion to Euler
sim_out.euler = zeros(length(sim_out.time), 3);
for i = 1:length(sim_out.time)
    sim_out.euler(i,:) = quatToEuler( sim_out.attitude(i,:)' );
end
sim_out.euler_d = zeros(length(sim_out.time), 3);
for i = 1:length(sim_out.time)
    sim_out.euler_d(i,:) = quatToEuler( sim_out.att_sp(i,:)' );
end

%clear i