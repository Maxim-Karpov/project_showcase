%Lac operon model - function
function dxdt=lac_rep_model(t,x)
 
%initialise model vector with zeros
dxdt=zeros(5,1);

%rate constant parameters
k1 = 10; %LacI + alac --> LacI-alac
k2 = 10; %LacI-alac --> LacI + alac
k3 = 1; %LacI + g_op --> LacI-g_op
k4 = 0.1; %LacI-g_op --> LacI + g_op

%ordinary differential equations for each reaction
dxdt(1) = k2*x(4)+k4*x(5)-k1*x(1)*x(2)-k3*x(1)*x(3); % d[LacI]/dt
dxdt(2) = k2*x(4)-k1*x(1)*x(2); % d[alac]/dt
dxdt(3) = k4*x(5)-k3*x(1)*x(3); % d[g_op]/dt
dxdt(4) = k1*x(1)*x(2)-k2*x(4); % d[LacI-alac]/dt
dxdt(5) = k3*x(1)*x(3)-k4*x(5); % d[LacI-g_op]/dt
