% Euler Method function
function [t_vals, y_vals] = euler_method(dydt, t0, y0, t_end, h)
    t_vals = t0:h:t_end;
    y_vals = zeros(size(t_vals));
    y_vals(1) = y0;
    for i = 1:length(t_vals)-1
        y_vals(i+1) = y_vals(i) + h * dydt(t_vals(i), y_vals(i));
    end
end

% % Taylor Method of order 2 function
% function [t_vals, y_vals] = taylor_method_order2(dydt, d2ydt2, t0, y0, t_end, h)
%     t_vals = t0:h:t_end;
%     y_vals = zeros(size(t_vals));
%     y_vals(1) = y0;
%     for i = 1:length(t_vals)-1
%         f = dydt(t_vals(i), y_vals(i));
%         f_prime = d2ydt2(t_vals(i), y_vals(i));
%         y_vals(i+1) = y_vals(i) + h * f + (h^2 / 2) * f_prime;
%     end
% end

% Modified Euler Method
function [t_vals, y_vals] = modified_euler(dydt, t0, y0, t_end, h)
    t_vals = t0:h:t_end;
    y_vals = zeros(size(t_vals));
    y_vals(1) = y0;
    for i = 1:length(t_vals) - 1
        y_mid = y_vals(i) + h / 2 * dydt(t_vals(i), y_vals(i));
        y_vals(i + 1) = y_vals(i) + h * dydt(t_vals(i) + h / 2, y_mid);
    end
end


% Taylor Method of order 2 function
function [t_vals, y_vals] = taylor_method_order2(dydt, d2ydt2, t0, y0, t_end, h)
    t_vals = t0:h:t_end;
    y_vals = zeros(size(t_vals));
    y_vals(1) = y0;
    for i = 1:length(t_vals)-1
        f = dydt(t_vals(i), y_vals(i));
        f_prime = d2ydt2(t_vals(i), y_vals(i));
        y_vals(i+1) = y_vals(i) + h * f + (h^2 / 2) * f_prime;
    end
end

% Taylor Method of order 4 function
function [t_vals, y_vals] = taylor_method_order4(dydt, d2ydt2, d3ydt3, d4ydt4, t0, y0, t_end, h)
    t_vals = t0:h:t_end;
    y_vals = zeros(size(t_vals));
    y_vals(1) = y0;
    for i = 1:length(t_vals)-1
        f = dydt(t_vals(i), y_vals(i));
        f_prime = d2ydt2(t_vals(i), y_vals(i));
        f_double_prime = d3ydt3(t_vals(i), y_vals(i));
        f_triple_prime = d4ydt4(t_vals(i), y_vals(i));
        y_vals(i+1) = y_vals(i) + h * f + (h^2 / 2) * f_prime + (h^3 / 6) * f_double_prime + (h^4 / 24) * f_triple_prime;
    end
end

% Runge-Kutta Method of order 4 function
function [t_vals, y_vals] = runge_kutta_order4(dydt, t0, y0, t_end, h)
    t_vals = t0:h:t_end;
    y_vals = zeros(size(t_vals));
    y_vals(1) = y0;
    for i = 1:length(t_vals)-1
        k1 = dydt(t_vals(i), y_vals(i));
        k2 = dydt(t_vals(i) + h/2, y_vals(i) + h/2 * k1);
        k3 = dydt(t_vals(i) + h/2, y_vals(i) + h/2 * k2);
        k4 = dydt(t_vals(i) + h, y_vals(i) + h * k3);
        y_vals(i+1) = y_vals(i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % QUESTION-1
% 
% % Define functions for each part of Question 1 and plot using Euler’s Method
% 
% % Question 1(a)
% dydt_1a = @(t, y) t * exp(3 * t) - 2 * y;
% exact_1a = @(t) (1/25) * (5 * t .* exp(3 * t) - exp(-2 * t)) + 1/25;
% [t_vals_1a, y_vals_1a] = euler_method(dydt_1a, 0, 0, 1, 0.5);
% 
% % Plot for 1(a)
% figure;
% subplot(2, 2, 1);
% hold on;
% title("Q1(a): Euler's Method");
% plot(t_vals_1a, y_vals_1a, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_1a, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 1(b)
% dydt_1b = @(t, y) 1 + (t - y)^2;
% exact_1b = @(t) t + 1/(1 - t);
% [t_vals_1b, y_vals_1b] = euler_method(dydt_1b, 2, 1, 3, 0.5);
% 
% % Plot for 1(b)
% subplot(2, 2, 2);
% hold on;
% title("Q1(b): Euler's Method");
% plot(t_vals_1b, y_vals_1b, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_1b, [2, 3], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 1(c)
% dydt_1c = @(t, y) 1 + y/t;
% exact_1c = @(t) t .* log(t) + 2 * t;
% [t_vals_1c, y_vals_1c] = euler_method(dydt_1c, 1, 2, 2, 0.25);
% 
% % Plot for 1(c)
% subplot(2, 2, 3);
% hold on;
% title("Q1(c): Euler's Method");
% plot(t_vals_1c, y_vals_1c, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_1c, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 1(d)
% dydt_1d = @(t, y) cos(2 * t) + sin(3 * t);
% exact_1d = @(t) (1/2) * sin(2 * t) - (1/3) * cos(3 * t) + 4/3;
% [t_vals_1d, y_vals_1d] = euler_method(dydt_1d, 0, 1, 1, 0.25);
% 
% % Plot for 1(d)
% subplot(2, 2, 4);
% hold on;
% title("Q1(d): Euler's Method");
% plot(t_vals_1d, y_vals_1d, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_1d, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %QUESTION-3
% 
% % Define functions for each part of Question 3 and plot using Euler’s Method
% 
% % Question 3(a)
% dydt_3a = @(t, y) y/t - (y/t)^2;
% exact_3a = @(t) t ./ (1 + log(t));
% [t_vals_3a, y_vals_3a] = euler_method(dydt_3a, 1, 1, 2, 0.1);
% 
% % Plot for 3(a)
% figure;
% subplot(2, 2, 1);
% hold on;
% title("Q3(a): Euler's Method");
% plot(t_vals_3a, y_vals_3a, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_3a, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 3(b)
% dydt_3b = @(t, y) 1 + y/t + (y/t)^2;
% exact_3b = @(t) t .* tan(log(t));
% [t_vals_3b, y_vals_3b] = euler_method(dydt_3b, 1, 0, 3, 0.2);
% 
% % Plot for 3(b)
% subplot(2, 2, 2);
% hold on;
% title("Q3(b): Euler's Method");
% plot(t_vals_3b, y_vals_3b, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_3b, [1, 3], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 3(c)
% dydt_3c = @(t, y) -(y + 1) * (y + 3);
% exact_3c = @(t) -3 + 2 ./ (1 + exp(-2 * t));
% [t_vals_3c, y_vals_3c] = euler_method(dydt_3c, 0, -2, 2, 0.2);
% 
% % Plot for 3(c)
% subplot(2, 2, 3);
% hold on;
% title("Q3(c): Euler's Method");
% plot(t_vals_3c, y_vals_3c, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_3c, [0, 2], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% % Question 3(d)
% dydt_3d = @(t, y) -5 * y + 5 * t^2 + 2 * t;
% exact_3d = @(t) t.^2 + (1/3) * exp(-5 * t);
% [t_vals_3d, y_vals_3d] = euler_method(dydt_3d, 0, 1/3, 1, 0.1);
% 
% % Plot for 3(d)
% subplot(2, 2, 4);
% hold on;
% title("Q3(d): Euler's Method");
% plot(t_vals_3d, y_vals_3d, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_3d, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %QUESTION-5
% 
% % Define function for Question 5(a) and apply Euler’s Method
% 
% dydt_5a = @(t, y) (2/t) * y + t^2 * exp(t);
% exact_5a = @(t) t.^2 .* (exp(t) - exp(1));
% [t_vals_5a, y_vals_5a] = euler_method(dydt_5a, 1, 0, 2, 0.1);
% 
% % Plot for 5(a)
% figure;
% hold on;
% title("Q5(a): Euler's Method");
% plot(t_vals_5a, y_vals_5a, 'ro-', 'DisplayName', 'Euler Approximation');
% fplot(exact_5a, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
% legend;
% hold off;
% 
% %part-(b)
% % Approximate values using linear interpolation
% y_approx_1_04 = interp1(t_vals_5a, y_vals_5a, 1.04);
% y_approx_1_55 = interp1(t_vals_5a, y_vals_5a, 1.55);
% y_approx_1_97 = interp1(t_vals_5a, y_vals_5a, 1.97);
% 
% % Display the interpolated results
% fprintf('Approximation of y(1.04): %.5f\n', y_approx_1_04);
% fprintf('Approximation of y(1.55): %.5f\n', y_approx_1_55);
% fprintf('Approximation of y(1.97): %.5f\n', y_approx_1_97);
% 
% % Compare with the exact values
% exact_y_1_04 = exact_5a(1.04);
% exact_y_1_55 = exact_5a(1.55);
% exact_y_1_97 = exact_5a(1.97);
% fprintf('Exact value of y(1.04): %.5f\n', exact_y_1_04);
% fprintf('Exact value of y(1.55): %.5f\n', exact_y_1_55);
% fprintf('Exact value of y(1.97): %.5f\n', exact_y_1_97);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUESTION-9

% Define function and second derivative for 9(a)
dydt_9a = @(t, y) t * exp(3 * t) - 2 * y;
d2ydt2_9a = @(t, y) exp(3*t) * (1 + t) + 4 * y;
exact_9a = @(t) (1/5) * t .* exp(3 * t) - (1/25) * exp(3 * t) + (1/25) * exp(-2 * t);

% Apply Taylor's Method of Order 2
[t_vals_9a, y_vals_9a] = taylor_method_order2(dydt_9a, d2ydt2_9a, 0, 0, 1, 0.5);

% Plot for 9(a)
figure;
subplot(2, 2, 1);
hold on;
title("Q9(a): Taylor's Method Order 2");
plot(t_vals_9a, y_vals_9a, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_9a, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 9(b)
dydt_9b = @(t, y) 1 + (t - y)^2;
d2ydt2_9b = @(t, y) -2 * (t - y)^3;
exact_9b = @(t) t + 1 ./ (1 - t);

% Apply Taylor's Method of Order 2
[t_vals_9b, y_vals_9b] = taylor_method_order2(dydt_9b, d2ydt2_9b, 2, 1, 3, 0.5);

% Plot for 9(b)
subplot(2, 2, 2);
hold on;
title("Q9(b): Taylor's Method Order 2");
plot(t_vals_9b, y_vals_9b, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_9b, [2, 3], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 9(c)
dydt_9c = @(t, y) 1 + y/t;
d2ydt2_9c = @(t, y) 1/t;
exact_9c = @(t) t .* log(t) + 2 * t;

% Apply Taylor's Method of Order 2
[t_vals_9c, y_vals_9c] = taylor_method_order2(dydt_9c, d2ydt2_9c, 1, 2, 2, 0.25);

% Plot for 9(c)
subplot(2, 2, 3);
hold on;
title("Q9(c): Taylor's Method Order 2");
plot(t_vals_9c, y_vals_9c, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_9c, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 9(d)
dydt_9d = @(t, y) cos(2 * t) + sin(3 * t);
d2ydt2_9d = @(t, y) -2 * sin(2 * t) + 3 * cos(3 * t);
exact_9d = @(t) (1/2) * sin(2 * t) - (1/3) * cos(3 * t) + (4/3);

% Apply Taylor's Method of Order 2
[t_vals_9d, y_vals_9d] = taylor_method_order2(dydt_9d, d2ydt2_9d, 0, 1, 1, 0.25);

% Plot for 9(d)
subplot(2, 2, 4);
hold on;
title("Q9(d): Taylor's Method Order 2");
plot(t_vals_9d, y_vals_9d, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_9d, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUESTION-11

% Define function and second derivative for 11(a)
dydt_11a = @(t, y) y/t - (y/t)^2;
d2ydt2_11a = @(t, y) -y^2/t^3 * (1 - 2*y/t);
exact_11a = @(t) t ./ (1 + log(t));

% Apply Taylor's Method of Order 2
[t_vals_11a, y_vals_11a] = taylor_method_order2(dydt_11a, d2ydt2_11a, 1, 1, 1.2, 0.1);

% Plot for 11(a)
figure;
subplot(2, 2, 1);
hold on;
title("Q11(a): Taylor's Method Order 2");
plot(t_vals_11a, y_vals_11a, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_11a, [1, 1.2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 11(b)
dydt_11b = @(t, y) sin(t) + exp(-t);
d2ydt2_11b = @(t, y) cos(t) - exp(-t);

% Apply Taylor's Method of Order 2
[t_vals_11b, y_vals_11b] = taylor_method_order2(dydt_11b, d2ydt2_11b, 0, 0, 1, 0.5);

% Plot for 11(b)
subplot(2, 2, 2);
hold on;
title("Q11(b): Taylor's Method Order 2");
plot(t_vals_11b, y_vals_11b, 'ro-', 'DisplayName', 'Taylor Approximation');

% Numerical approximation of exact solution via integral
integral_vals = arrayfun(@(t) integral(@(tau) sin(tau) + exp(-tau), 0, t), t_vals_11b);
plot(t_vals_11b, integral_vals, 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 11(c)
dydt_11c = @(t, y) (y^2 + y) / t;
d2ydt2_11c = @(t, y) 2*y/t * (y^2/t + y/t);
exact_11c = @(t) -3 + 2 ./ (1 + exp(-2 * t));

% Apply Taylor's Method of Order 2
[t_vals_11c, y_vals_11c] = taylor_method_order2(dydt_11c, d2ydt2_11c, 1, -2, 3, 0.5);

% Plot for 11(c)
subplot(2, 2, 3);
hold on;
title("Q11(c): Taylor's Method Order 2");
plot(t_vals_11c, y_vals_11c, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_11c, [1, 3], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define function and second derivative for 11(d)
dydt_11d = @(t, y) -t * y + 4 * t / y;
d2ydt2_11d = @(t, y) -y - t * dydt_11d(t, y) +4/y - 4*t * dydt_11d(t, y)/y^2;
exact_11d = @(t) t.^2 + (1/3) * exp(-5 * t);

% Apply Taylor's Method of Order 2
[t_vals_11d, y_vals_11d] = taylor_method_order2(dydt_11d, d2ydt2_11d, 0, 1, 1, 0.25);

% Plot for 11(d)
subplot(2, 2, 4);
hold on;
title("Q11(d): Taylor's Method Order 2");
plot(t_vals_11d, y_vals_11d, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_11d, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUESTION-13

% Define first and second derivatives for 13(a)
dydt_13 = @(t, y) (2 / t) * y + t^2 * exp(t);
d2ydt2_13 = @(t, y) (2 / t) * dydt_13(t, y) + 2 * exp(t) + t^2 * exp(t);

% Exact solution
exact_13 = @(t) t.^2 .* (exp(t) - exp(1));

% Apply Taylor's Method of Order 2
[t_vals_13a, y_vals_13a] = taylor_method_order2(dydt_13, d2ydt2_13, 1, 0, 2, 0.1);

% Plot for 13(a)
figure;
subplot(2, 1, 1);
hold on;
title("Q13(a): Taylor's Method Order 2");
plot(t_vals_13a, y_vals_13a, 'ro-', 'DisplayName', 'Taylor Approximation');
fplot(exact_13, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Linear interpolation for values in 13(b)
y_approx_1_04 = interp1(t_vals_13a, y_vals_13a, 1.04);
y_approx_1_55 = interp1(t_vals_13a, y_vals_13a, 1.55);
y_approx_1_97 = interp1(t_vals_13a, y_vals_13a, 1.97);

% Display results
fprintf('Approximation of y(1.04) from 13(a): %.5f\n', y_approx_1_04);
fprintf('Approximation of y(1.55) from 13(a): %.5f\n', y_approx_1_55);
fprintf('Approximation of y(1.97) from 13(a): %.5f\n', y_approx_1_97);

% Exact values for comparison
exact_y_1_04 = exact_13(1.04);
exact_y_1_55 = exact_13(1.55);
exact_y_1_97 = exact_13(1.97);
fprintf('Exact y(1.04): %.5f\n', exact_y_1_04);
fprintf('Exact y(1.55): %.5f\n', exact_y_1_55);
fprintf('Exact y(1.97): %.5f\n', exact_y_1_97);

% Define third and fourth derivatives for Taylor's method of order 4
d3ydt3_13 = @(t, y) (2/t) * d2ydt2_13(t, y) - 4 * dydt_13(t, y) / t^2 + (2 * exp(t) + t^2 * exp(t));
d4ydt4_13 = @(t, y) (2/t) * d3ydt3_13(t, y) - 6 * d2ydt2_13(t, y) / t^2 + (2 * exp(t) + t^2 * exp(t));

% Apply Taylor's Method of Order 4
[t_vals_13c, y_vals_13c] = taylor_method_order4(dydt_13, d2ydt2_13, d3ydt3_13, d4ydt4_13, 1, 0, 2, 0.1);

% Plot for 13(c)
subplot(2, 1, 2);
hold on;
title("Q13(c): Taylor's Method Order 4");
plot(t_vals_13c, y_vals_13c, 'bo-', 'DisplayName', 'Taylor Approximation (Order 4)');
fplot(exact_13, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Linear interpolation for values in 13(d)
y_approx_1_04_4 = interp1(t_vals_13c, y_vals_13c, 1.04);
y_approx_1_55_4 = interp1(t_vals_13c, y_vals_13c, 1.55);
y_approx_1_97_4 = interp1(t_vals_13c, y_vals_13c, 1.97);

% Display results
fprintf('Approximation of y(1.04) from 13(c): %.5f\n', y_approx_1_04_4);
fprintf('Approximation of y(1.55) from 13(c): %.5f\n', y_approx_1_55_4);
fprintf('Approximation of y(1.97) from 13(c): %.5f\n', y_approx_1_97_4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUESTION-15

% Define the function for 15(a)
dydt_15a = @(t, y) y/t - (y/t)^2;
exact_15a = @(t) t ./ (1 + log(t));

% Apply Modified Euler and Runge-Kutta methods
[t_vals_15a_me, y_vals_15a_me] = modified_euler(dydt_15a, 1, 1, 2, 0.1);
[t_vals_15a_rk, y_vals_15a_rk] = runge_kutta_order4(dydt_15a, 1, 1, 2, 0.1);

% Plot for 15(a)
figure;
subplot(2, 2, 1);
hold on;
title("Q15(a): Modified Euler and RK4");
plot(t_vals_15a_me, y_vals_15a_me, 'ro-', 'DisplayName', 'Modified Euler');
plot(t_vals_15a_rk, y_vals_15a_rk, 'bo-', 'DisplayName', 'Runge-Kutta (Order 4)');
fplot(exact_15a, [1, 2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define the function for 15(b)
dydt_15b = @(t, y) 1 + y/t + (y/t)^2;
exact_15b = @(t) t .* tan(log(t));

% Apply Modified Euler and Runge-Kutta methods
[t_vals_15b_me, y_vals_15b_me] = modified_euler(dydt_15b, 1, 0, 3, 0.2);
[t_vals_15b_rk, y_vals_15b_rk] = runge_kutta_order4(dydt_15b, 1, 0, 3, 0.2);

% Plot for 15(b)
subplot(2, 2, 2);
hold on;
title("Q15(b): Modified Euler and RK4");
plot(t_vals_15b_me, y_vals_15b_me, 'ro-', 'DisplayName', 'Modified Euler');
plot(t_vals_15b_rk, y_vals_15b_rk, 'bo-', 'DisplayName', 'Runge-Kutta (Order 4)');
fplot(exact_15b, [1, 3], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define the function for 15(c)
dydt_15c = @(t, y) -(y + 1) * (y + 3);
exact_15c = @(t) -3 + 2 ./ (1 + exp(-2 * t));

% Apply Modified Euler and Runge-Kutta methods
[t_vals_15c_me, y_vals_15c_me] = modified_euler(dydt_15c, 0, -2, 2, 0.2);
[t_vals_15c_rk, y_vals_15c_rk] = runge_kutta_order4(dydt_15c, 0, -2, 2, 0.2);

% Plot for 15(c)
subplot(2, 2, 3);
hold on;
title("Q15(c): Modified Euler and RK4");
plot(t_vals_15c_me, y_vals_15c_me, 'ro-', 'DisplayName', 'Modified Euler');
plot(t_vals_15c_rk, y_vals_15c_rk, 'bo-', 'DisplayName', 'Runge-Kutta (Order 4)');
fplot(exact_15c, [0, 2], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;

% Define the function for 15(d)
dydt_15d = @(t, y) -5 * y + 5 * t^2 + 2 * t;
exact_15d = @(t) t.^2 + (1/3) * exp(-5 * t);

% Apply Modified Euler and Runge-Kutta methods
[t_vals_15d_me, y_vals_15d_me] = modified_euler(dydt_15d, 0, 1/3, 1, 0.1);
[t_vals_15d_rk, y_vals_15d_rk] = runge_kutta_order4(dydt_15d, 0, 1/3, 1, 0.1);

% Plot for 15(d)
subplot(2, 2, 4);
hold on;
title("Q15(d): Modified Euler and RK4");
plot(t_vals_15d_me, y_vals_15d_me, 'ro-', 'DisplayName', 'Modified Euler');
plot(t_vals_15d_rk, y_vals_15d_rk, 'bo-', 'DisplayName', 'Runge-Kutta (Order 4)');
fplot(exact_15d, [0, 1], 'k--', 'DisplayName', 'Exact Solution');
legend;
hold off;