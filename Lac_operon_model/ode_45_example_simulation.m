#example simulation
%ode_45
% Settings
options = odeset('InitialStep',0.1,'MaxStep',0.1);
t_range= [0 5];
x_ini= [15 0 10 0 0];
 
% Simulation
[t,x]=ode45(@lac_rep_model,t_range,x_ini,options);
 
% Plot time series of all variables
figure(1);
plot(t,x(:,:));
xlabel('Time');
ylabel('Species concentration');
legend('LacI','alac','g_o_p','LacI-alac','LacI-g_o_p','Location','NorthEast');
