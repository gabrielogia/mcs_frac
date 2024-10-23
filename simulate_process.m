function W = simulate_process(time, ndof, EPS, nonstat, is_base, frequency_max)

    dt = time(2) - time(1);
    N = numel(time);
    %frequency_max = pi/dt;
    
    m=800;
    dw = frequency_max/m;
        
    W = zeros(N,ndof); 

    if is_base

        if nonstat
            for ii=0:m-1
                %EPS = S0*((ii*dw/5/pi)^2).*exp(-0.1.*(time)).*(time.^2).*exp(-((ii*dw/5/pi)^2).*time);
                eps_fun = EPS{1}.fun(ii*dw,time);
                W(:,1) = W(:,1) + (((4*dw.*eps_fun).^0.5).*cos((ii*dw).*time - (2*pi)*rand))'; 
            end
        else
            S0 = EPS{1}.fun(0,0);
            W(:,1) = sqrt(2*pi*S0/dt).*randn(N,1);
        end

        for jj=2:ndof
    
            W(:,jj) = W(:,1);
    
        end

    else

        for jj=1:ndof
    
            if nonstat
                for ii=0:m-1
                    %EPS = S0*((ii*dw/5/pi)^2).*exp(-0.1.*(time)).*(time.^2).*exp(-((ii*dw/5/pi)^2).*time);
                    eps_fun = EPS{jj}.fun(ii*dw,time);
                    W(:,jj) = W(:,jj) + (((4*dw.*eps_fun).^0.5).*cos((ii*dw).*time - (2*pi)*rand))'; 
                end
            else
                S0 = EPS{jj}.fun(0,0);
                W(:,jj) = sqrt(2*pi*S0/dt).*randn(N,1);
            end
    
    
        end

    end

    W = W';

end
     