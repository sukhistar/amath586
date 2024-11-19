% HW2_C1
% Sukhjit Kaur
% AMATH586

clc
close
clear

e = exp(1); %Euler's number e
x_0 = -5; %x_0
x_end = 5; %final x
Nx = [20, 40, 80, 160];%Nx to try
T = 3;
t0 = 0; %start time
tf = T; %final time


%Explicit Euler (Forward)

for i = 1:length(Nx)
    N = Nx(i);
    dx = (x_end-x_0)/N;
    dt = 0.4*(dx)^2;
    t = t0:dt:tf;


    %initial condition
    % u(t,x) t => m, x => j convention 
    % u(0,x) = v(1,x)

    for j = 1:N
        for m = 1:length(t)
            v(m,j) = 1/(sqrt(4*pi*m))*e^(-(j-2)^2/4*m);
        end
        u(1,j) = v(2,j);%initial condition  u(0,x) = v(1,x)
    end

    %exact solution

     for j=1:N
               u_exact(1,j) = v(1,j);
         u_exact(end,j) = v(end,j);
    for m = 1:length(t)-1
  
            u_exact(m,j) = v(m+1,j);
        end
     end
   

    %boundary conditions
    % u(t,-5) = v(t+1,-5)
    % u(t,5) = v(t+1,5)
    for m = 2:length(t)-1
        u(m,1) = v(m+1,1);
        u(m, N) = v(m+1,N);
    end
    %implement forward euler
    lambda = dt/(dx)^2;
    for m = 1:length(t)
        for j = 2:N-1
            u(m+1,j) = lambda*(u(m,j+1)-2*u(m,j)+u(m,j-1))+u(m,j);

        end
    end
    FE = u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Implicit Euler (Backward)

for i = 1:length(Nx)
    N = Nx(i);
    dx = (x_end-x_0)/N;
    dt = dx;
    t = t0:dt:tf;

    %initial condition
    % u(t,x)
    % u(0,x) = v(1,x)

    for j = 1:N
        for m = 1:length(t)-1
            v(m,j) = 1/(sqrt(4*pi*m))*e^(-(j-2)^2/4*m);
        end
        u(1,j) = v(2,j);%initial condition  u(0,x) = v(1,x)
    end

    %boundary conditions
    % u(t,-5) = v(t+1,-5)
    % u(t,5) = v(t+1,5)
    for m = 1:length(t)-2
        u(m,1) = v(m+1,1);
        u(m, N) = v(m+1,N);
    end
    %implement backward euler
    lambda = dt/(dx)^2;
    for m = 1:length(t)-1
        for j = 2:N-1
            u(m+1,j) = lambda*(u(m+1,j+1)-2*u(m+1,j)+u(m+1,j-1))+u(m,j);
        end
    end
    BE = u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Crank-Nicolson

for i = 1:length(Nx)
    N = Nx(i);
    dx = (x_end-x_0)/N;
    dt = dx;
    t = t0:dt:tf;

    %initial condition
    % u(t,x)
    % u(0,x) = v(1,x)

    for j = 1:N
        for m = 1:length(t)-1
            v(m,j) = 1/(sqrt(4*pi*m))*e^(-(j-2)^2/4*m);
        end
        u(1,j) = v(2,j);%initial condition  u(0,x) = v(1,x)
    end

    %boundary conditions
    % u(t,-5) = v(t+1,-5)
    % u(t,5) = v(t+1,5)
    for m = 1:length(t)-2
        u(m,1) = v(m+1,1);
        u(m, N) = v(m+1,N);
    end
    %implement crank-nicolson scheme
    lambda = dt/2*(dx)^2;
    for m = 1:length(t)-1
        for j = 2:N-1
            u(m+1,j) = lambda*(u(m,j+1)-2*u(m,j)+u(m,j-1)+u(m+1,j+1)-2*u(m+1,j)+u(m+1,j-1))+u(m,j);
        end
    end
    CN = u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate errors
%Explicit Euler - expect O(dt+dx^2) = O(dx^2)
FE_error = norm(u_exact-FE)/norm(u_exact);
%Implicit Euler - expect O(dt+dx^2) = O(dx)
BE_error = norm(u_exact-BE)/norm(u_exact);
%Crank-Nicolson - expect O(dt^2+dx^2) = O(dx^2)
CN_error = norm(u_exact-CN)/norm(u_exact);