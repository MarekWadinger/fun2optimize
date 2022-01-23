function [x, xlog, iter] = gdm(f, x0, kmax, bounds, opts)
    arguments
        f
        x0 (1,:) double {mustBeNumeric, mustBeReal} = ''
        kmax (1,1) int16 {mustBeNumeric, mustBeReal} = 500
        bounds (1,2) double {mustBeNumeric, mustBeReal} = [0, 10]
        opts.t (1,1) double {mustBeNumeric, mustBeReal} = 0.2
        opts.Plot (1,1) logical = false
        opts.Printout (1,1) logical = true
    end
    tic
    
    if opts.Plot clf(figure(101)); figure(101); end

    syms a
    df = str2func(['@(a)' char(diff(f(a)))]);
    
    if x0
        x = x0;
    else
        % Define Neighborhood
        d = bounds(2) - bounds(1);
        % Pick initial value from the neighborhood
        x = bounds(1) + rand()*d;
    end

    xlog = zeros(length(x), kmax+1);
    xlog(:,1) = x;
    
    for iter = 1:kmax
        % Compute negative of the gradient in x
        dx = -df(x);
    
        % Stopping criterion - gradient sufficiently close to 0
        if norm(dx) <= 10^(-6)
            break
        end
        
        % Update optimized variable
        x = x+opts.t*dx;
        xlog(:, iter+1) = x;
    
        if opts.Plot
            subplot(121); plot(iter, x, 'b.'); title('x');
            drawnow limitrate; hold on
            subplot(122); plot(iter, f(x), 'b.'); title('f(x)');
            drawnow limitrate; hold on
        end
    
    end
    
    if opts.Printout
        xlog = xlog(:,1:iter);
        text = strcat('Optimum reached in: ', num2str(iter), ...
            ' iterations... Time: ', num2str(toc), 's');
        disp(text)
    end
end
