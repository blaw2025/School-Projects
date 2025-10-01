clear; clc; close all;

%% Full Procedure

%  Collect boundary points
num_points = 1000;  % 1000 sample points in x
x_space = linspace(-2, 1, num_points);
y_space = nan(size(x_space)); % Initialize Nan array

for k = 1:num_points
    fn = indicator_fn_at_x(x_space(k));

    try
		y_space(k) = bisection(fn, 0, 2);% lower bound s=0, upper bound e=2
    catch
        y_space(k) = NaN;
    end
end

% Plot raw boundary samples
figure;
plot(x_space, y_space, 'b.');
xlabel('Re(c)'); 
ylabel('Im(c)');
title('Mandelbrot boundary samples');
axis equal tight;

% Fit polynomial
valid = ~isnan(y_space); % Discard NaN and flat tails

% Trim x-range to avoid flat tails 
in_range = x_space >= -2 & x_space <= 0.5; 

points = valid & in_range; % Keep desired points
x_fit = x_space(points); % x points to fit
y_fit = y_space(points); % y points to fit

% Fit polynomial of order 15 with polyfit
order = 15;
p = polyfit(x_fit, y_fit, order);

% Evaluate polynomial at new x values
x = linspace(min(x_fit), max(x_fit)); % Create x range
y = polyval(p, x);

% Plot fit
figure;
plot(x_space, y_space, 'b.'); 
hold on;
plot(x, y, 'r-', 'LineWidth', 1.5);
xlabel('Re(c)'); 
ylabel('Im(c)');
title(sprintf('Order %d Polynomial', order));
legend('Boundary samples','Polynomial fit','Location','best');
axis equal tight;

% Compute Boundary length
s = min(x_fit); % Left bound on x
e = max(x_fit); % Right bound on x
L = poly_len(p, s, e); % Curve length

% Print 2*L for full length
fprintf('Approximate total circumference ≈ %.6f\n', 2*L);

%% function it = fractal(c)

% Computes number of iterations till divergence
function it = fractal(c)
    z = 0;
    for k = 1:100 % Choose 100 max iterations
        if abs(z) > 2
            it = k;  % Number of iterations
            return; % End function and return the number of iterations
		end
		z = z^2 + c; % Update z accordingly
    end
    it = 0; % Return 0 if c doesn't diverge
end

%% function m = bisection(fn_f, s, e)

% Finds the point on the boundary of the fractal.

% Indicator function 
function fn = indicator_fn_at_x(x)
    fn = @(y) (fractal(x + 1i*y) > 0)*2 - 1; % Focus on a line
end

% Bisection method
function m = bisection(fn_f, s, e)
    fs = fn_f(s); % y-value of fn_f with x-value s
    fe = fn_f(e); % y-value of fn_f with x-value e

	% Print error if there's no sign change
    if fs * fe > 0
        error('No sign change in [s,e]');
    end
	
	% Otherwise, perform bisection method
    for k = 1:100 % Choose 100 max iterations
        m = 0.5*(s+e); % Start from the middle of s and e
        fm = fn_f(m);

		% End method if tolerance is met
        if (e - s)/2 < 1e-6 % Stopping condition tolerance
            return;
		end

		% Checks which side m is on
        if fs * fm <= 0 % Negative side
            e = m; % Update lower bound and start from a new middle 
		else % Otherwise, positive side
            s = m; % Update upper bound and start from a new middle 
			fs = fm;
        end
    end
    m = 0.5*(s+e); % Return the average of s and e
end

%% function l = poly_len(p, s, e)

% Compute curve length of a polynomial

function l = poly_len(p, s, e)
    dp = polyder(p); % Matlab's derivative function
    dydx = @(x) polyval(dp, x); % Make dy/dx a math function
    ds = @(x) sqrt(1 + (dydx(x)).^2); % Arc length to be integrated
    l = integral(ds, s, e); % Matlab's integral function
end
