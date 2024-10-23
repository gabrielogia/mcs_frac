clc
close all

ns=1;
epsilon=ex;
number_of_samples = ns;
    % non-linear parameter
    %e = epsilon;

    % fractional
    q = alpha;

    Mi = inv(M);
    MiC = Mi*C;
    MiK = Mi*K;
    
    % Excitation
    p1 = 5*pi; 
    p2 = 5*pi;
    for i=1:ndof
        EPS{i}.fun = @(f,t)( S0*((f./p1).^2).*exp(-b0.*t).*(t.^2).*exp(-((f./p2).^2).*t));
    end

    % Nonstationary: true
    nonstat = true;
    
    % Number of samples
    ns = number_of_samples;
    
    %MCS
    tt = time;
    
    for i=1:ns
        w = simulate_process(tt, ndof, EPS, nonstat, is_base);
    
        xini = zeros(3*ndof,1);
        [t1,X] = fde45_mdof(@(t,x) FUNODE_mdof(t,x,MiC,MiK,Mi,stiffness,epsilon,ndof,w,tt), [tt(1) tt(end)], xini, q);

        state(i, :, :) = X(1:ndof, :); %just the displacement
        stored_time(i, :) = t1;
    end

    time_out = stored_time(1,:);

    [~, ~, n_time] = size(state); % z time dimension.

    for i=1:n_time
        for j=1:ndof
            vx(i, j) = var(state(:, j, i));
        end
    end

    %vx = interp1(t1, vx, time);

    vx = vx';


    %% ==========================
    % First passage time

    first_passage_time = zeros(ns,ndof);
    for i=1:ns
        
        for j=1:ndof
            sample_path = state(i,j,:);
            
            time_aux = time_out(abs(sample_path) > barrier);
            first_passage_time(i,j) = time_aux(1);
        end

    end
        

