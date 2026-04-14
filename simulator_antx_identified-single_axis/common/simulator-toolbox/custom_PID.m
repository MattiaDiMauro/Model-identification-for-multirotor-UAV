    function [ ctr ] = custom_PID( Kff, Kp, Ki, Kd, N, ts )

feedForward_term = tf(Kff);
feedForward_term.u = 'y_0';
feedForward_term.y = 'u_ff';

proportional_term = tf(Kp);
proportional_term.u = 'e';
proportional_term.y = 'u_p';

integrator_tf = tf(ts, [1, -1], ts);
integral_term = Ki * integrator_tf;
integral_term.u = 'e';
integral_term.y = 'u_i';

%Derivative action with CASCADED low pass filter
% derivative_tf = tf([1, -1], [ts, 0], ts);
% derivative_term = Kd * derivative_tf;
% derivative_term.u = 'y';
% derivative_term.y = 'u_d';
%
% low_pass_filter = tf([ts, 0], [(1/N + ts), -1/N], ts); %Backward Euler
% % low_pass_filter = tf(ts, [1/N, (ts - 1/N)], ts); %Forward Euler
% low_pass_filter.u = 'u_d';
% low_pass_filter.y = 'u_d_f';

%Derivative action with BUILT-IN low pass filter
high_pass_filter_tf = tf([1, -1], [(1/N + ts), -1/N], ts);
derivative_term = Kd * high_pass_filter_tf;
derivative_term.u = 'y';
derivative_term.y = 'u_d_f';

%% Blocks sum
error_sum = sumblk('e = y_0 - y');
control_variable_sum = sumblk('u = u_ff + u_p + u_i - u_d_f');

%% Controllers
ctr = connect(feedForward_term, ...
    proportional_term, integral_term, derivative_term, ...
    error_sum, control_variable_sum, ...
    {'y_0', 'y'},{'e', 'u'});

end