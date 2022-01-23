function [x, xlog, iter] = luusjaakola(f, x0, kmax, bounds, opts)
% Luus-Jaacola algorithm
    arguments
        f
        x0 double {mustBeNumeric, mustBeReal} = ''
        kmax (1,1) int16 {mustBeNumeric, mustBeReal} = 500
        bounds (1,2) double {mustBeNumeric, mustBeReal} = [0, 10]
        opts.ShrinkRate (1,1) double {mustBeNumeric, mustBeReal, ...
            mustBeInRange(opts.ShrinkRate,0,1)} = 0.98
        opts.Plot (1,1) logical = false
        opts.Printout (1,1) logical = true
    end
    tic
    
    if opts.Plot clf(figure(101)); figure(101); end
    
    % Define Neighborhood
    d = bounds(2) - bounds(1);

    if x0;
        x = x0;
    else
        % Pick initial value from the neighborhood
        x = bounds(1) + rand()*d;
    end

    xlog = [];

    for iter = 1:kmax
        % Pick random vector
        a = -d/2 + rand()*d; 
        % Set new candidate
        y = x + a;
        
        % Boundaries represent whole domain of the function
        while y >=bounds(2) | y<bounds(1)
            a = -d/2 + rand()*d;
            y = x + a;
        end
    
        % Stopping criterion - change of obj fcn sufficiently small
        if abs(f(y) - f(x))<1e-6
            break
        end
        
        % Update candidate if obj fcn is smaller
        if f(y) < f(x)
            x = y;
            xlog = [xlog; y];
        % Otherwise shrink neighborhood
        else
            d = opts.ShrinkRate * d;
        end
    
        if opts.Plot
            subplot(121); plot(iter, x, 'b.'); title('x');
            drawnow limitrate; hold on
            subplot(122); plot(iter, f(x), 'b.'); title('f(x)');
            drawnow limitrate; hold on
        end
    
    end
    
    if opts.Printout
        text = strcat('Optimum reached in: ', num2str(iter), ...
            ' iterations... Time: ', num2str(toc), 's');
        disp(text)
    end
end
