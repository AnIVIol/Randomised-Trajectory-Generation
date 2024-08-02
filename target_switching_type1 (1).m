%% Basics
clc;
clear all;
close all;

%% Global constants
global v targets_x targets_y f r_max;
v = 1;
time = [0, 400];
r_min = 3;
r_max = 5;
r_range = linspace(r_min, r_max, 100);
max_range = linspace(0, r_max);
targets_x = [-r_min, r_min, r_min, -r_min];
targets_y = [r_min, r_min, -r_min, -r_min];

%% Generating g(r)
coeff = gen_gr(v, r_min, r_max);
a = coeff(1); b = coeff(2); c = coeff(3); d = coeff(4);
g = @(r) a.*(r.^3) + b.*(r.^2) + c.*r + d;
f = @(r) 3*a*r + 2*b + c./r;

%% Generating trajectory
tol = 1e-13;
options = odeset('RelTol', tol, 'AbsTol', tol);
global prev_target current_target
prev_target=2;
current_target= 1 ;
[t1, Y1] = ode45(@(t1,Y1) odefunc(t1, Y1), time, [targets_x(1)-r_max, targets_y(1), -pi/2], options);

%% Plotting
figure(1);
        hold on;
        grid on;
    plot(max_range, -v*max_range, 'k');
    plot(max_range, v*max_range, 'k');
    plot(r_range, g(r_range), 'b');
        hold off;

angle = 0:0.01:2*pi;
figure(2);
        set(gcf, 'Position', [1000 100 800 800]);
        hold on;
        grid on;
    for i = 1:length(targets_x)
        plot(targets_x(i) + r_min*cos(angle), targets_y(i) + r_min*sin(angle), '--b');
        plot(targets_x(i) + r_max*cos(angle), targets_y(i) + r_max*sin(angle), '--b');
    end
    plot(targets_x, targets_y, 'xk');
    comet(Y1(:,1), Y1(:,2))
        axis equal;
        hold off;

%% ODE function
function ret = odefunc(t, Y)
    global v f r_max;
    global targets_x targets_y;
    global prev_target
    global current_target

    % Unpacking variables
    x = Y(1);
    y = Y(2);
    alpha = Y(3);

    % Unicycle dynamics
    x_dot = v * cos(alpha);
    y_dot = v * sin(alpha);

    r_temp = r_max;
    current_temp = current_target;
    for i = 1:length(targets_x)
        
        x_2 = (targets_x(i) - x)^2;
        y_2 = (targets_y(i) - y)^2;
        if(sqrt(x_2+y_2)<r_temp)
            if(i~=prev_target)
                current_target = i;
                r_temp = sqrt(x_2+y_2);
            end
        end
        
    end
    r_ = r_temp;
    if(current_target~=current_temp)
        prev_target = current_temp;
    end
    alpha_dot = f(r_);

    ret = [x_dot; y_dot; alpha_dot];
end

%% Functions used to generate generating functions
function ret = gen_gr(v, r_min, r_max)
    n = rand();
    
    % Linear function generation
    if (n < 1)
        c = v*(r_max+r_min) / (r_max-r_min);
        d = (v-c) * r_max;
        ret = [0, 0, c, d];
    end
end