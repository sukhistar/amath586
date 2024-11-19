% AMATH 586
% HW 3 - HR w/ mc
% Sukhjit Kaur

clear;
clc;

%Variables
a = 1;
x_min = 0;
x_max = 2;
N = 200;
cfl = 0.8;
t = 0;
t_final = 1;

theta = zeros(1,N+3);
F_rl = zeros(1,N+3);
F_rh =zeros(1,N+3);
F_ll = zeros(1,N+3);
F_lh = zeros(1,N+3);
F_right = zeros(1,N+3);
F_left  = zeros(1,N+3);

%Discretize
dx = (x_max-x_min)/N;
x = x_min - dx : dx : x_max + dx;
dt = cfl*dx;
tspan = t_final/dt;

%Initial condition
u0 = zeros(1,N+3);
w = find ( 0 <= x & x <= 0.6 );
u0(w) = exp(-100*(x(w)-0.3).^2);
j = find (0.8 <= x & x<= 1);
u0(j) = 1;
u = u0;
unext = zeros(1,N+3);
ue = u0;

% solve
for k = 1:tspan
    for i = 2:N-1

if u(i) == u(i+1)
    theta(i) = 1;
else
    theta(i) = (u(i) - u(i-1)) / (u(i+1) - u(i)); 
    end
        theta(1) = 1; 
        theta(N) = 1;
    end

 phi = max(0, min(min( (1+theta)/2, 2), 2*theta)); %mc
 for i = 2:N+2    
        % Compute fluxes
        F_rl(i) = u(i);
        F_rh(i) = (1/2)*(u(i)+u(i+1)) - (1/2)*(cfl)*(u(i+1)-u(i));
        
        F_ll(i) = u(i-1);
        F_lh(i) = (1/2)*(u(i-1)+u(i)) - (1/2)*(cfl)*(u(i)-u(i-1));
        
        % Compute next time step
        F_right(i) = F_rl(i) + phi(i)*( F_rh(i) - F_rl(i) );
        F_left(i)  = F_ll(i) + phi(i-1)*( F_lh(i) - F_ll(i) );
        
        unext(i) = u(i) - cfl*(F_right(i) - F_left(i));
 end
        u = unext;
        t=t+dt;
        roundt = round(t,2);

    % Calculate exact solution
    ue = zeros(1,N+3);
    w = find ( 0 <= (x-t) & (x-t) <= 0.6 );
    ue(w) = exp(-100*(x(w)-0.3-t).^2);
    j = find (0.8 <= (x-t) & (x-t)<= 1);
    ue(j) = 1;

    % %Plot at t=0.5
    if roundt == 0.50
    plot(x,ue,'b-', 'markerfacecolor', 'b')
    hold on
    plot(x,u,'r-', 'MarkerFaceColor','r')
    hold off
    axis([x_min x_max -0.5 1.5])
    xlabel('x','FontSize',16)
    ylabel('U(t,x)','FontSize',16)
    title(sprintf('time = %1.3f',t))
    legend('Exact Solution', 'High Resolution with mc limiter')
    shg
    pause(0.5)
    end
     % %Plot at t=0.8
    if roundt == 0.80
        figure
    plot(x,ue,'b-', 'markerfacecolor', 'b')
    hold on
    plot(x,u,'r-', 'MarkerFaceColor','r')
    hold off
    axis([x_min x_max -0.5 1.5])
    xlabel('x','FontSize',16)
    ylabel('U(t,x)','FontSize',16)
    title(sprintf('time = %1.3f',t))
    legend('Exact Solution', 'High Resolution with mc limiter')
    shg
    pause(0.5)
    end
end
%  plot(x,u)
% hold on 
% plot(x,ue)
% hold off

