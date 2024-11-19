% AMATH 586
% HW 3 - HR w/ mc
% Sukhjit Kaur

clear;
clc;

%Variables
x_min = 0;
x_max = 2;
N = 200;
cfl = 0.8;
t = 0;
t_final = 1;

%Discretize
dx = (x_max-x_min)/N;
x = x_min - dx : dx : x_max + dx;
dt = cfl*dx;
tspan = t_final/dt;

%Initial conditions
u0 = zeros(1,N+3);
w = find ( 0 <= x & x <= 0.6 );
u0(w) = exp(-100*(x(w)-0.3).^2);
j = find (0.8 <= x & x<= 1);
u0(j) = 1;
u = u0;
ue = u0;

% solve
for k = 1:tspan
    for i = 2:N

        if u(i) ~= 0
            r(i) = (u(i) - u(i-1))/((u(i+1) - u(i)));
        end

        phi = max([0, min([2*r,0.5*(r+1),2])]); %mc

        Fl = u(i) + 0.5*(1-cfl)*(u(i+1)-u(i)).*phi;
        Fr = circshift(Fl,1);
        unext = u - cfl*(Fl - Fr);
    end
   
    u = unext;
    t = t + dt;
    roundt = round(t,2);

    % Calculate exact solution
    ue = zeros(1,N+3);
    w = find ( 0 <= (x-t) & (x-t) <= 0.6 );
    ue(w) = exp(-100*(x(w)-0.3-t).^2);
    j = find (0.8 <= (x-t) & (x-t)<= 1);
    ue(j) = 1;

    %Plot at t=0.5
    if roundt == 0.50
    plot(x,ue,'b-')
    hold on
    plot(x,u,'r-*', 'MarkerFaceColor','r')
    hold off
    axis([x_min x_max -0.5 1.5])
    xlabel('x','FontSize',16)
    ylabel('U(t,x)','FontSize',16)
    title(sprintf('time = %1.3f',t))
    legend('Exact Solution', 'High Resolution with mc limiter')
    shg
    pause(0.5)
    end

    %Plot at t=0.8
    if roundt == 0.80
        figure
    plot(x,ue,'b-')
    hold on
    plot(x,u,'r-*', 'MarkerFaceColor','r')
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