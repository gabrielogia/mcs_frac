function amp = get_amplitude(x)

    [r,c] = size(x);

    if r>c
        x = x';
    end

    npadd = ceil(max(r,c)/10);
    % X = padarray(x',npadd,'circ')';
    X = padarray(x',npadd,'symmetric')';

    if r>c
        X = X';
    end

  
    A = abs(hilbert(X));
    amp = A(npadd+1:end-npadd);


    
    
