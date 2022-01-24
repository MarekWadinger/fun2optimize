function [x, fval, xlog, iter] = gdmb(f, x0, kmax, bounds, opts)
    arguments
        f
        x0 double {mustBeNumeric, mustBeReal} = 0
        kmax (1,1) int16 {mustBeNumeric, mustBeReal} = 500
        bounds (:,2) double {mustBeNumeric, mustBeReal} = [0, 10]
        opts.t (1,1) double {mustBeNumeric, mustBeReal} = 1
        opts.alpha (1,1) double {mustBeNumeric, mustBeReal} = 0.2
        opts.beta (1,1) double {mustBeNumeric, mustBeReal} = 0.5
        opts.Plot (1,1) logical = false
        opts.Printout (1,1) logical = true
    end
    tic
    
    if opts.Plot clf(figure(101)); figure(101); end

    if x0
        % assign number of optimized variables
        nopts = length(x0);

        x = x0;
    else
        % assign number of optimized variables
        nopts = size(bounds,1);
        % Define Neighborhood
        d = bounds(:,2) - bounds(:,1);
        % Pick initial value from the neighborhood
        x = (bounds(:,1) + rand(nopts,1).*d)';
    end
    
    % create symbolic variable
    syms xsym [1 nopts]
    % create character vector of symbolic variables
    charopts = char(xsym);
    % if more than one, remove square brackets
    if nopts > 1 
        charopts = charopts(2:end-1);
    end
    
    % create cell array
    xsymcell = sym2cell(xsym);
    
    % Create function handle for gradient and hessian
    df = str2func(['@(' charopts ')' ...
                  char( gradient( f(xsymcell{:}) ) ) ]);

    % Allocate memory for x logging
    xlog = zeros(kmax+1, length(x));
    xlog(1,:) = x;

    % Assign optimization parameters
    t = opts.t; alpha = opts.alpha; beta = opts.beta;

    for iter = 1:kmax
        % Compute negative of the gradient in x
        xcell = num2cell(x);
        dx = -df(xcell{:});
        
        % Stopping criterion - gradient sufficiently close to 0
        if norm(dx) <= 10^(-6)
            break
        end
    
        xnewcell = num2cell(x+opts.t*dx');
        % Backtracking - update step size
        while f(xnewcell{:}) > f(xcell{:})+alpha*t*df(xcell{:})'*dx
            t = beta*t;
            xnewcell = num2cell(x+t*dx');
        end

        % Update optimized variable
        x = x+t*dx';
        xlog(iter+1, :) = x;
    
        if opts.Plot
            subplot(121); plot(iter, x, 'b.'); title('x');
            drawnow limitrate; hold on
            subplot(122); plot(iter, f(x), 'b.'); title('f(x)');
            drawnow limitrate; hold on
        end
    
    end
    
    xlog = xlog(1:iter, :);

    if opts.Printout
        text = strcat('Optimum reached in: ', num2str(iter), ...
            ' iterations... Time: ', num2str(toc), 's');
        disp(text)
    end

    fval = f(xnewcell{:});
end
