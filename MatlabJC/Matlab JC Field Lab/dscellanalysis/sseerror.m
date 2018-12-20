function sse = sseerror(func, x, y, params)
% data: column 1: contrast; column 2: response
    sse = 0;
    for i = 1:length(x)
        pred = func(params,x(i));
        error = (pred-y(i)).^2;
        sse = error + sse;
    end
end