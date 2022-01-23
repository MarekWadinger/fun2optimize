function [x, xlog, iter] = simulatedannealing(f, x0, kmax, bounds, opts)
    arguments
        f
        x0 double {mustBeNumeric, mustBeReal} = 0
        kmax (1,1) double {mustBeNumeric, mustBeReal} = 500
        bounds (1,2) double {mustBeNumeric, mustBeReal} = [0, 10]
        opts.Temp (1,1) double {mustBeNumeric, mustBeReal} = 50
        opts.T_reduction (1,1) double {mustBeNumeric, mustBeReal, ...
            mustBeInRange(opts.T_reduction,0,1)} = 0.98
        opts.Plot (1,1) logical = false
        opts.Printout (1,1) logical = true
    end
    tic
    
    if opts.Plot clf(figure(102)); figure(102); end
    
    % Define Neighborhood
    d = bounds(2) - bounds(1);
    
    if x0
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
        
        % Compute Boltzmann Probability
        P = min([1 exp(-(f(y)-f(x))/opts.Temp)]);

        % Update candidate if Probability is 1 - value of obj fcn is smaller
        % or value of Boltzmann P is bigger than random number from uniform
        % distribution.
        if P == 1 | P > rand
            x = y;                  % assign to x
            xlog = [xlog; y];
        end
        
        % Decrease temperature
        opts.Temp = opts.Temp * opts.T_reduction;
    
        if opts.Plot
            subplot(221); plot(iter, 1, 'bo');
            xlim([0,kmax]); title('terminal iteration');
            drawnow limitrate
            subplot(222); plot(iter, opts.Temp, 'b.'); title('Temperature');
            drawnow limitrate; hold on
            subplot(223); plot(iter, x, 'b.'); title('x');
            drawnow limitrate; hold on
            subplot(224); plot(iter, f(x), 'b.'); title('f(x)');
            drawnow limitrate; hold on
        end
    
    end
    
    if opts.Printout
        text = strcat('Optimum reached in: ', num2str(iter), ...
            ' iterations... Time: ', num2str(toc), 's');
        disp(text)
    end
end