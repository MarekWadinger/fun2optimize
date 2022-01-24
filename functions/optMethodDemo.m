function [x, f_val, t, i] = optMethodDemo(f, x0, kmax, bounds, ...
    opt_method, opt_param, opt_values)
    arguments
        f
        x0 double {mustBeNumeric, mustBeReal} = 0
        kmax (1,1) double {mustBeNumeric, mustBeReal} = 500
        bounds (:,2) double {mustBeNumeric, mustBeReal} = [0, 10]
        opt_method string = 'gdm'
        opt_param string = 't'
        opt_values (1,:) double {mustBeNumeric, mustBeReal} = 1:1:10
    end

    x = zeros(length(opt_values),1);
    f_val = zeros(length(opt_values),1);
    t = zeros(length(opt_values),1);
    i = zeros(length(opt_values),1);
    
    opt_method_fcn = str2func(opt_method)

    for t_curr = 1:length(opt_values)
        tic
        [xn, ~, iter] = opt_method_fcn(f, x0, kmax, bounds, ...
                                   opt_param, opt_values(t_curr) , ...
                                   'Plot', false, ...
                                   'Printout', false);
        x(t_curr) = xn;
        f_val(t_curr) = -f(xn);
        t(t_curr) = toc;
        i(t_curr) = iter;
    end
    
    clf(figure(3)); figure(3)
    subplot(221); plot(opt_values, x, 'b*'); title('Optimized variable value'); 
    subplot(222); plot(opt_values, f_val, 'b*'); title('Objective function value');
    subplot(223); plot(opt_values, t, 'b*'); title('Time');
    subplot(224); plot(opt_values, i, 'b*'); title('Iterations');

    explainParam(opt_method, opt_param)
end