clc
clear
close all

rng(121)

% The parameters are defined in the following code:

r2 = 0.05; % intensity of parametric excitation
r = sqrt(r2);
Ca = 0.1; % Fractional damping coefficient.
b0=Ca;
w02 = 10; % Natural frequency squared.
w0 = sqrt(w02); % Natural frequency.

S0=b0*w02/pi; % White noise magnitude.

q = 0.5; % Fractional derivative.
epx = 0.5; % Nonlinearity parameter.

ndof = 2;
M = [1 0;1 1];
C = [0.2 -0.2;
    0 0.2];
K = [10 -10;
    0 10];

Mi = inv(M);
MiC = Mi*C;
MiK = Mi*K;
gfun = @(x)(x);

% Power Spectrum
%EPS{1}.fun = @(f,t)( 0.1*((f./p1).^2).*exp(-0.05.*t).*(t.^2).*exp(-((f./p2).^2).*t));nonstat = true;
nonstat = false;
for i=1:ndof
    EPS{i}.fun = @(f,t)( S0.*ones(size(f)).*ones(size(t)));
end


%MCS
tf = 5;
dt = 0.001;
tt = 0:dt:tf;
nt = numel(tt);
ntlim = nt-floor(nt/15);

% Number of MC simulations.
ns = 1;

fprintf('Number of simulations: %d\n', ns);
fprintf('Order fractional derivative: %f\n', q);
fprintf('Nonlinearity parameter: %d\n', epx);

tic
for i=1:ns
    w = simulate_process(tt, ndof, EPS,nonstat);

    xini = zeros(3*ndof,1);
    % Solve the Fractional differential equation.
    % [t1,X] = fde45_mdof(@(t,x) FUNODE_mdof(t,x,b0,w02,epx,r,w,u,tt), [tt(1) tt(end)], xini, q);
    [t1,X] = fde45_mdof(@(t,x) FUNODE_mdof(t,x,MiC,MiK,Mi,gfun,ndof,w,tt), [tt(1) tt(end)], xini, q);
  
end
elapsed_time = toc;

Filename = sprintf('Data_Ex3_%s.mat', datestr(now,'mm-dd-yyyy HH-MM'));

fprintf('Elapsed time: %f\n', elapsed_time);
fprintf('File name: %s\n', Filename);

%save(Filename)


