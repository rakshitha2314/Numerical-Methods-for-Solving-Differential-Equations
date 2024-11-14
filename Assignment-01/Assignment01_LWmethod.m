% MTL 712 

% ASSIGNMENT 01 - Part 01 - LF

% Objective: For the given 5 test cases, plot the graphs using Lax -
% Freirdrich method.




% fromula to be used:
% u(i, j+1) = 0.5 * (u(i+1, j) + u(i-1, j)) - lambda/2 * (u(i+1, j) -
% u(i-1, j))



function [] = lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)
    delta_x = (x_final - x_init) / x_parts;
    delta_t = lambda * delta_x;
    t_parts = (t_req - t_init)/delta_t;
    x = linspace(x_init, x_final, x_parts);
    
    u = init_cond;
    f_u = zeros(x_parts);
    for i=1:x_parts
        f_u(i) = f(u(i));
    end

    %time-stepping using Lax-Friedrichs:
    for n = 1:t_parts
        u_new = zeros(size(init_cond));
        for i = 2: x_parts - 1
            u_new(i) = u(i) - 0.5 * lambda * (f_u(i+1) - f_u(i-1)) + 0.5 * lambda^2 * (f_u(i+1) - 2 * f_u(i) + f_u(i-1));
        end
        % Periodic boundary conditions:
        u_new(1) = u(1) - 0.5 * lambda * (f_u(2) - f_u(x_parts)) + 0.5 * lambda^2 * (f_u(2) - 2 * f_u(1) + f_u(x_parts));
        u_new(x_parts) = u(x_parts) - 0.5 * lambda * (f_u(1) - f_u(x_parts - 1)) + 0.5 * lambda^2 * (f_u(1) - 2 * f_u(x_parts) + f_u(x_parts - 1));
        % Update variable 
        u = u_new;
        f_u = zeros(x_parts);
        for i=1:x_parts
            f_u(i) = f(u(i));
        end
    end

    % Generate required plots:
    figure;
    plot(x, exact_soln, 'k', 'LineWidth', 2); % Exact Solution
    hold on;
    plot(x, u, 'ko', 'LineWidth', 2, 'MarkerSize', 5); % Numerical Solution
    title(sprintf('Lax-Wendroff Method at t = %d', t_req));
    xlabel('x-val');
    ylabel('u-req');
    % legend('Exact Solution', 'Lax-Friedrichs Solution');
    grid on;
end


function [] = lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)
    delta_x = (x_final - x_init) / x_parts;
    delta_t = lambda * delta_x;
    t_parts = (t_req - t_init)/delta_t;
    x = linspace(x_init, x_final, x_parts);
    
    u = init_cond;
    f_u = zeros(x_parts);
    for i=1:x_parts
        f_u(i) = f(u(i));
    end

    %time-stepping using Lax-Friedrichs:
    for n = 1:t_parts
        u_new = zeros(size(init_cond));
        for i = 2: x_parts - 1
            u_new(i) = 0.5 * (u(i-1) + u(i+1)) - 0.5 * lambda * (f_u(i+1) - f_u(i-1)) + 0.5 * lambda^2 * (f_u(i+1) - 2 * f_u(i) + f_u(i-1));
        end
        % Periodic boundary conditions:
        u_new(1) = 0.5*(u(2) + u(x_parts)) - 0.5 * lambda * (f_u(2) - f_u(x_parts)) + 0.5 * lambda^2 * (f_u(2) - 2 * f_u(1) + f_u(x_parts));
        u_new(x_parts) = 0.5*(u(x_parts-1) + u(1)) - 0.5 * lambda * (f_u(1) - f_u(x_parts - 1)) + 0.5 * lambda^2 * (f_u(1) - 2 * f_u(x_parts) + f_u(x_parts - 1));
        % Update variable 
        u = u_new;
        f_u = zeros(x_parts);
        for i=1:x_parts
            f_u(i) = f(u(i));
        end
    end

    % Generate required plots:
    figure;
    plot(x, exact_soln, 'k', 'LineWidth', 2); % Exact Solution
    hold on;
    plot(x, u, 'ko', 'LineWidth', 2, 'MarkerSize', 5); % Numerical Solution
    title(sprintf('Lax-Wendroff Method at t = %d', t_req));
    xlabel('x-val');
    ylabel('u-req');
    % legend('Exact Solution', 'Lax-Friedrichs Solution');
    grid on;
end

% % Test Case - 01
% 
% x_init = -1;
% x_final = 1; 
% x_parts = 40; 
% lambda = 0.8; 
% t_init = 0;
% t_req = 30;
% 
% x = linspace(x_init, x_final, x_parts);
% init_cond = -sin(pi*x); 
% 
% exact_soln = init_cond;
% f = @(x) x;
% 
% lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% 
% 
% % Test Case - 02
% 
% x_init = -1;
% x_final = 1; 
% x_parts = 40; 
% lambda = 0.8; 
% t_init = 0;
% t_req = 4;
% 
% x = linspace(x_init, x_final, x_parts);
% init_cond = zeros(x_parts);
% for i = 1:x_parts
%     if abs(x(i)) < 1/3
%         init_cond(i) = 1;
%     else
%         init_cond(i) = 0;
%     end
% end
% 
% exact_soln = init_cond;
% f = @(x) x;
% 
% lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% 
% 
% 
% 
% % Test Case - 03
% 
% x_init = -1;
% x_final = 1; 
% x_parts = 600; 
% lambda = 0.8; 
% t_init = 0;
% t_req_1 = 4;
% t_req_2 = 40;
% 
% x = linspace(x_init, x_final, x_parts);
% init_cond = zeros(x_parts);
% for i = 1:x_parts
%     if abs(x(i)) < 1/3
%         init_cond(i) = 1;
%     else
%         init_cond(i) = 0;
%     end
% end
% 
% exact_soln = init_cond;
% f = @(x) x;
% 
% lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req_1, init_cond, exact_soln,f)
% lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req_1, init_cond, exact_soln,f)
% 
% 
% lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req_2, init_cond, exact_soln,f)
% lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req_2, init_cond, exact_soln,f)





% Test Case - 04

x_init = -1;
x_final = 1; 
x_parts = 40; 
lambda = 0.8; 
t_init = 0;
t_req = 0.6;


x = linspace(x_init, x_final, x_parts);
init_cond = zeros(x_parts);
for i = 1:x_parts
    if abs(x(i)) > 1/2
        init_cond(i) = 0;
    else
        init_cond(i) = -2*abs(x(i)) + 1;
    end
end

exact_soln = init_cond;
% for i = 1:x_parts
%     if x(i) > -1/3 && x(i) < -1/3 + t_req
%         exact_soln(i) = (x(i) + 1/3)/t_req;
%     elseif x(i) < 1/3 + t_req/2 && x(i) > -1/3 + t_req
%         exact_soln(i) = 1;
%     end
% end
        

f = @(x) x^2/2;

lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
%lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)



% % Test Case - 05
% 
% x_init = -1;
% x_final = 1; 
% x_parts = 40; 
% lambda = 0.8; 
% t_init = 0;
% t_req = 0.3;
% 
% 
% x = linspace(x_init, x_final, x_parts);
% init_cond = zeros(x_parts);
% for i = 1:x_parts
%     if abs(x(i)) < 1/3
%         init_cond(i) = 1;
%     else
%         init_cond(i) = -1;
%     end
% end
% 
% exact_soln = zeros(x_parts);
% for i = 1:x_parts
%     if x(i) < -1/3 - t_req
%         exact_soln(i) = -1;
%     elseif x(i) > -1/3 - t_req && x(i) < -1/3 + t_req
%         exact_soln(i) = -1 + (x(i) + 1/3 + t_req)/t_req;
%     elseif x(i) < 1/3 && x(i) > -1/3 + t_req
%         exact_soln(i) = 1;
%     else
%         exact_soln(i) = -1;
%     end
% end
% 
% 
% f = @(x) x^2/2;
% 
% lax_wendroff(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% lax_wendroff_modified_for_average(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)
% 
% 
