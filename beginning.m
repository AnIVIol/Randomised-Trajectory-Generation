 clc;
clear all;
close all;

global v;
v = 1;
time = [0, 100];
r_min = 2;
r_max = 3;
K = -15;
r_range = linspace(r_min, r_max, 100);
max_range = linspace(0, r_max);

m = v * (r_max+r_min) / (r_max-r_min);
c = v*r_max - m*r_max;

syms a  b;
eqn1 = @(a,b)  a*(r_min)^2+b*r_min == +v*r_min-K;
eqn2 = @(a,b) a*(r_max)^2+b*r_max == -v*r_max-K;
sol = solve({eqn1, eqn2}, [a, b]);

global g f;
%g = @(r) m.*r + c;
%f = @(r) -m./r;
g = @(r) (sol.a).*r.^2+(sol.b).*r+K;
f = @(r) 2*(sol.a)+(sol.b)./r;

tol = 1e-5;
options = odeset('RelTol', tol, 'AbsTol', tol);
[t1, Y1] = ode45(@(t1,Y1) odefunc(t1, Y1), time, [-r_max, 0, pi/2], options);

%% Plotting
figure(2);
        hold on;
        grid on;
    plot(r_range, g(r_range), 'b');
    plot(max_range, v*max_range, 'k');
    plot(max_range, -v*max_range, 'k');
        hold off;

angle = 0:0.01:2*pi;
figure(1);
        set(gcf, 'Position', [1000 100 800 800]);
        hold on;
        grid on;
    plot(r_min*cos(angle), r_min*sin(angle), '--b');
    plot(r_max*cos(angle), r_max*sin(angle), '--b');
%     plot(Y1(:,1), Y1(:,2), 'b')
    comet(Y1(:,1), Y1(:,2));
        hold off;


%% Functions used
function ret = odefunc(t, Y)
    global f v;

    % Unpacking variables
    x = Y(1);
    y = Y(2);
    alpha = Y(3);
    r_ = sqrt(x^2 + y^2);

    % Unicycle dynamics
    x_dot = v * cos(alpha);
    y_dot = v * sin(alpha);
    alpha_dot = f(r_);

    % Return
    ret = [x_dot; y_dot; alpha_dot];
end