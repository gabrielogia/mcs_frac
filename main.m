%% structure data
clc
clear
close all
%%
% For reproducibility:
rng(1111);

% Is Base motion / non-stationary (excitation):
is_base = false;
nonstat = true;

% Number of DOFs:
ndof = 2;

% Fractional derivative:
q = 0.5; 

% Nonlinearity parameter:
epx = 0.2*ones(1,ndof);

% Mass, damping, and stiffness vectors: 
mass = 1*ones(1,ndof); 
damping = 50*ones(1,ndof);
stiffness = 100*ones(1,ndof);

% Maximum time:
T = 10;

% Barrier:
lam = 1;

% Time increment for the Monte Carlo simulation.
dT = 0.001; %dT = 0.0001;

% Construct matrices M, C, and K:
[M, C, K] = get_mck(mass, damping, stiffness, ndof);

% Maximum frequency of the power spectrum:
fmax_ps = 50; 

% Number of samples in the MCS:
ns = 48;

% Discretization in time and frequency for the Statistical Linearization:
ntime = 200;
nfreq = 1000;

% Run MCS:
run_mcs = true;

% Find the Survival Probability (Yes=1, No=0):
run_fps = true;

%% Statistical Linearization
disp("Running Statistical Linearization:")

time = linspace(0, T, ntime);
omega_n = sqrt(eig(inv(M)*K));
freq = linspace(0,fmax_ps,nfreq);

[varx_sl, varv_sl, conv, k_eq, c_eq] = ...
    statistical_linearization(mass, damping, stiffness, M, C, K,freq, time, ndof, epx, q, is_base, oscillator);

for i=1:ndof
    smaxt(i) = max(varx_sl(1,:));
end

smaxi = max(smaxt);
for i=1:ndof
    
    barrier(i) = lam*sqrt(smaxi);
end

%% Monte Carlo Simulation
%run_mcs = input('Run Monte Carlo Simulation (Yes=1, No=0):');

if run_mcs
    disp('Running MCS:')
    
    [varx_mcs, time_out, first_passage_time,response,amplitude] = ...
        monte_carlo(ns,M,C,K,epx,q,mass,damping,stiffness,fmax_ps, nonstat, is_base,T,dT, barrier);
end

%%

idd=1;
figure
cc=1;
for i=1:ns
    z = amplitude(idd,:,i);
    %z = response(idd,:,i);
    B = barrier(idd);
    tb=time_out(abs(z)>B);
    if ~isempty(tb)
    xx=time_out;
    z(time_out>tb(1))=[];
    xx(time_out>tb(1))=[];
    hold on
    plot(xx,z,'k')
    plot(xx(end),B,'r.')
    tf(cc)=xx(end);
    cc=cc+1;
    end
end


%% Equivalent damping and stiffness
omega_eq_2 = varv_sl./varx_sl;

for i=1:ndof

    omega_eq_2(i,1) = findfirstpoint(omega_eq_2(i,2),omega_eq_2(i,3));

    for j=1:numel(time)
        t=time(j);
        sig2t = varx_sl(i,j);
        Sw = @(x)( evolutionary_power_spectrum(x, t) );
        %Sx = @(x,y)( Sw(x)./( (omega_eq_2(i,j) - x.^2).^2 + (y*x).^2 )  );
        Sx = @(x,y)( Sw(x)./( abs(omega_eq_2(i,j) - x.^2 + y*(1i*x).^q).^2 )  );
        %sfun = @(y) ( (sig2t - pi*omega_eq_2(i,j)./( omega_eq_2(i,j)*y ) ).^2  );

        sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );

        bt = fminbnd(sfun,0.01,1500);
        beq(i,j)=bt;

    end
    %beq(i,1) = beq(i,2);
    beq(i,1) = findfirstpoint(beq(i,2),beq(i,3));

end

beta_eq = beq;

%beta_eq = beta_eq./(sqrt(omega_eq_2).^(q-1).*sin(q*pi/2));
%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")


ic = 0.00000001;
c = zeros(ndof, numel(time));
for i=1:ndof
    beta_eq_dof = beta_eq(i,:);
    omega_eq_2_dof = omega_eq_2(i,:);

    beta_eq_dof(1) = findfirstpoint(beta_eq_dof(2),beta_eq_dof(3));
    omega_eq_2_dof(1) = findfirstpoint(omega_eq_2_dof(2),omega_eq_2_dof(3));

    [t, c(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q), time, ic);
end


%% Get c(t) by solving the ODE from stochastic averaging.
disp("Solving the ODE to find c(t):")

nitera=0;
%figure
cn=1;
if nitera>0
    c = zeros(ndof, numel(time));
end

niterav = 1:ndof;
for i=1:ndof
    for ii=1:nitera

        ic = 0.00000001;
        beta_eq_dof = beta_eq(i,:);
        omega_eq_2_dof = omega_eq_2(i,:);
    
        beta_eq_dof(1) = findfirstpoint(beta_eq_dof(2),beta_eq_dof(3));
        omega_eq_2_dof(1) = findfirstpoint(omega_eq_2_dof(2),omega_eq_2_dof(3));
    
        [t, c(i,:)] = ode89(@(t, c_aux) solve_c_mdof(t, c_aux, beta_eq_dof, omega_eq_2_dof, time,q), time, ic);
        ints2(i,ii) = trapz(time,c(i,:));

        for j=1:numel(time)
            t=time(j);
            sig2t = c(i,j);
            Sw = @(x)( evolutionary_power_spectrum(x, t) );
            %Sx = @(x,y)( Sw(x)./( (omega_eq_2(i,j) - x.^2).^2 + (y*x).^2 )  );
            Sx = @(x,y)( Sw(x)./( abs(omega_eq_2(i,j) - x.^2 + y*(1i*x).^q).^2 )  );
            sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );
    
            bt = fminbnd(sfun,0.000000001,1500);
            beq(i,j)=bt;
    
        end

        beq(i,1) = findfirstpoint(beq(i,2),beq(i,3));
        beta_eq(i,:) = beq(i,:);
    end
end
% plot(time,dc','r')
% hold on
% drawnow


%%
% j=3;
% i=1;
% t=time(j);
% sig2t = c(i,j);
% Sw = @(x)( evolutionary_power_spectrum(x, t) );
% Sx = @(x,y)( Sw(x)./( (omega_eq_2(i,j) - x.^2).^2 + (y*x).^2 )  );
% sfun = @(y) ( (sig2t - 2*integral(@(x)Sx(x,y),0,Inf)).^2  );
% 
% bt = fminbnd(sfun,0.000000001,1500);
% beq(i,j)=bt;
% clear x
% Y=0.0001:0.01:500;
% for i=1:numel(Y) 
%     x(i)=sfun(Y(i));
% end
% 
% close all
% semilogy(Y,x)

%%
figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time,omega_eq_2(i,:),'linewidth',2) % MCS
    xlabel('Time','interpreter','latex')
    ylabel('$\omega^2_{eq}(t)$','interpreter','latex')
end

figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, beta_eq(i,:),'linewidth',2)
    xlabel('Time','interpreter','latex')
    ylabel('$\beta_{eq}(t)$','interpreter','latex')
end

%%
figure('color',[1 1 1]);
for i=1:ndof
    subplot(ndof,1,i); 
    hold on
    plot(time, varx_sl(i,:),'k-','linewidth',2) % SL
    plot(time, c(i,:)','r--','linewidth',2) % ODE
    if run_mcs
        plot(time_out,varx_mcs(i,:),'b:','linewidth',2) % MCS
    end
    legend('SL', 'SA', 'MCS','interpreter','latex')
    xlabel('Time','interpreter','latex')
    ylabel('$Var[x(t)]$','interpreter','latex')
end

%% Survival Probability

%run_fps = input('Find the Survival Probability (Yes=1, No=0):');

if run_fps
    %bt = beta_eq./(sqrt(omega_eq_2).^(q-1).*sin(q*pi/2));
    bt = beta_eq;
    P=survival_probability(barrier,c,time,1000,bt,10);

    
    figure('color',[1 1 1]);
    for i=1:ndof
        fpt = first_passage_time(:,i);
        fpt = fpt(fpt>0);
        [fpp,tfp]=ksdensity(fpt,'width',0.1,'Function','survivor');
        subplot(ndof,1,i); 
        hold on
        scatter(fpt,linspace(0,1,numel(fpt)),10,'r','filled','MarkerFaceAlpha',0.2);
        plot(time, P(i,:)','k','linewidth',2);
        plot(tfp, fpp,'b--','linewidth',2);
        
        legend('MCS','Analytical','KDE','interpreter','latex')
        title('First passage time')
        xlabel('Time')
        ylabel('Survivor Propability')
        xlim([0 T])
        ylim([0 1])
    end

end