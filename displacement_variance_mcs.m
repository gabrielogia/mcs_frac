function vx = displacement_variance_mcs(omega0, epsilon, beta0, alpha, S0, b0, tmax, number_of_samples)
    w02 = omega0.^2; 
    epx = epsilon;
    Ca = beta0;
    
    % fractional
    q = alpha;
    
    % Excitation
    p1 = 5*pi; 
    p2 = 5*pi;
    EPS{1}.fun = @(f,t)( S0*((f./p1).^2).*exp(-b0.*t).*(t.^2).*exp(-((f./p2).^2).*t));
    
    % Nonstationary: true
    nonstat = true;
    
    % Number of samples
    ns = number_of_samples;
    
    % Final time instant
    tf = tmax;
    
    % Time discretization
    tt = 0:0.01:tf;
    tinter = 0:0.1:tf;
    
    nt = numel(tinter);
    Amplitude = zeros(ns,nt);
    X = zeros(ns,nt);
    
    parfor i=1:ns
        disp(i)
    
        % Simulate the excitation.
        w = simulate_process(tt, 1, EPS,nonstat);
    
        % Function handle to the ODE function.
        fun = @(t,x) FUNODE(t,x,w02,Ca,epx,w,tt);
    
        % initial conditions.
        xini = [0;0;0]; % 2*num_dof + num_fractional_derivatives (=num_dof)
    
        % Solve the ODE.
        [t1f,Xf] = fde45(fun, [tt(1) tt(end)], xini, q);
    
        % interpolate the displacemen
        x1 = interp1(t1f,Xf(:,1),tinter,'pchip');
    
        % Get amplitude (to plot the full PDF)
        Amplitude(i,:) = get_amplitude(x1);
    
        % Displacement.
        X(i,:) = x1;
    
    end
    
    % Get variance.
    vx = get_variance(X);

