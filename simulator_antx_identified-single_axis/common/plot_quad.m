%% Plot OUTPUT
figure
subplot(2,2,1)
hold on
plot(sim_out.time, sim_out.pos_ned);
set(gca,'ColorOrderIndex',1)
plot(sim_out.time, sim_out.pos_sp,'--');
grid minor
xlabel('[s]')
ylabel('[m]')
title('Position')
legend('N', 'E', 'D')

subplot(2,2,3)
plot(sim_out.time, sim_out.pos_sp - sim_out.pos_ned);
grid minor
xlabel('[s]')
ylabel('[m]')
title('Position error')
legend('N_{error}', 'E_{error}', 'D_{error}')

subplot(2,2,2)
hold on
plot(sim_out.time, sim_out.vel_ned);
set(gca,'ColorOrderIndex',1)
plot(sim_out.time, sim_out.vel_sp,'--');
grid minor
xlabel('[s]')
ylabel('[m/s]')
title('Velocity')
legend('N_{dot}', 'E_{dot}', 'D_{dot}')

subplot(2,2,4)
plot(sim_out.time, sim_out.vel_sp - sim_out.vel_ned);
grid minor
xlabel('[s]')
ylabel('[m/s]')
title('Velocity error')
legend('N_{dot error}', 'E_{dot error}', 'D_{dot error}')
 
%%
figure
subplot(2,2,1)
hold on
plot(sim_out.time, rad2deg(sim_out.euler));
set(gca,'ColorOrderIndex',1)
plot(sim_out.time, rad2deg(sim_out.euler_d),'--');
grid minor
xlabel('[s]')
ylabel('[deg]')
title('Attitude')
legend('R', 'P', 'Y')

subplot(2,2,3)
plot(sim_out.time, rad2deg(unwrap(sim_out.euler_d) - unwrap(sim_out.euler)));
grid minor
xlabel('[s]')
ylabel('[deg]')
title('Attitude error')
legend('R_{error}', 'P_{error}', 'Y_{error}')

subplot(2,2,2)
hold on
plot(sim_out.time, rad2deg(sim_out.ang_rates_body));
set(gca,'ColorOrderIndex',1)
plot(sim_out.time, rad2deg(sim_out.ang_vel_sp),'--');
grid minor
xlabel('[s]')
ylabel('[deg/s]')
title('Angular speed')
legend('p', 'q', 'r') 

subplot(2,2,4)
plot(sim_out.time, rad2deg(sim_out.ang_vel_sp - sim_out.ang_rates_body));
grid minor
xlabel('[s]')
ylabel('[deg/s]')
title('Angular speed error')
legend('p_{error}', 'q_{error}', 'r_{error}') 

%%
figure
plot(sim_out.time, sim_out.motors);
grid minor
xlabel('[s]')
ylabel('[%]')
title('Th% motors')
legend('Th% 1','Th% 2','Th% 3','Th% 4')

%% END OF CODE