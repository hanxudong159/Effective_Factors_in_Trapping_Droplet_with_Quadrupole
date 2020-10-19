clear;close;clc
Varray = [2000, 3000, 4000];
warray = [900, 1200, 1900];
for i = 1:3
    V = Varray(1);                                                      %magnitude of AC voltage
    w = warray(i);                                                      %frequency of AC voltage
    trange = [0, 0.5];                                                  %time range
    x0 = [0.001, 0.001];                                                %initial x position and velocity
    y0 = [0.001, 0.001];                                                %initial y position and velocity
    [t1, x] = ode45(@(t, x) xsimulation(t, x, V, w), trange, x0);   	%solve ode equation on x direction
    [t2, y] = ode45(@(t, y) ysimulation(t, y, V, w), trange, y0);   	%solve ode equation on y direction
    figure;                                                             %new figure
    plot(t1, x(:, 1));                                                  %plot x-t
    xlabel('t(s)');
    ylabel('x(m)');
    figure;                                                             %new figure
    plot(t2, y(:, 1));                                                  %plot y-t
    xlabel('t(s)');
    ylabel('y(m)');
    figure;                                                             %new figure
    t = linspace(0, 0.5, size(t1, 1) * 2);                              %new t linspace for interpolation
    xvals = interp1(t1, x(:, 1), t);                                    %x table after interpolation
    yvals = interp1(t2, y(:, 1), t);                                    %y table after intrepolation
    plot(xvals, yvals);                                                 %plot y-x for path of droplet after simulation
    xlabel('x(m)');
    ylabel('y(m)');
end

for i = 1:3
    V = Varray(i);                                                      %magnitude of AC voltage
    w = warray(2);                                                      %frequency of AC voltage
    trange = [0, 0.5];                                                  %time range
    x0 = [0.001, 0.001];                                                %initial x position and velocity
    y0 = [0.001, 0.001];                                                %initial y position and velocity
    [t1, x] = ode45(@(t, x) xsimulation(t, x, V, w), trange, x0);       %solve ode equation on x direction
    [t2, y] = ode45(@(t, y) ysimulation(t, y, V, w), trange, y0);       %solve ode equation on y direction
    figure;                                                             %new figure
    plot(t1, x(:, 1));                                                  %plot x-t
    xlabel('t(s)');
    ylabel('x(m)');
    figure;                                                             %new figure
    plot(t2, y(:, 1));                                                  %plot y-t
    xlabel('t(s)');
    ylabel('y(m)');
    figure;                                                             %new figure
    t = linspace(0, 0.5, size(t1, 1) * 2);                              %new t linspace for interpolation
    xvals = interp1(t1, x(:, 1), t);                                    %x table after interpolation
    yvals = interp1(t2, y(:, 1), t);                                    %y table after intrepolation
    plot(xvals, yvals);                                                 %plot y-x for path of droplet after simulation
    xlabel('x(m)');
    ylabel('y(m)');
end

function rk = xsimulation(t, x, V, w)
d = 50 * 10^(-6);                                                       %diameter of droplet
n = 1.849 * 10^(-5);                                                    %air dynamic viscosity
m = 4 / 3 * pi * (d / 2)^3 * (0.9974456 * 10^3);                        %mass of droplet
surface_tension = 7.28 * 10^(-2);                                       %surface tension
q = 0.3 * 8 * pi * sqrt(8.85e-12 * surface_tension * (d / 2)^3);        %maximum charge by Rayleigh limit
r0 = 1.2 * 10^(-2);                                                     %distance from center to rod surface
E = -2 * x(1) / r0^2 * V * cos(w * t);                                  %electric field on x
rk = zeros(2, 1);
rk(1) = x(2);                                                           %1st derivative of x
rk(2) = (E * q - 3 * pi * n * d * x(2)) / m;                            %2nd derivative of x
end

function rk = ysimulation(t, y, V, w)
d = 50 * 10^(-6);                                                       %diameter of droplet
n = 1.849 * 10^(-5);                                                    %air dynamic viscosity
m = 4 / 3 * pi * (d / 2)^3 * (0.9974456 * 10^3);                        %mass of droplet
surface_tension = 7.28 * 10^(-2);                                       %surface tension
q = 0.3 * 8 * pi * sqrt(8.85e-12 * surface_tension * (d / 2)^3);        %maximum charge by Rayleigh limit
r0 = 1.2 * 10^(-2);                                                     %distance from center to rod surface
E = 2 * y(1) / r0^2 * V * cos(w * t);                                   %electric field on y
rk = zeros(2, 1);
rk(1) = y(2);                                                           %1st derivative of y
rk(2) = (E * q - 3 * pi * n * d * y(2)) / m;                            %2nd derivative of y
end
