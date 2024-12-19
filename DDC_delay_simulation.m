%%paper: Data-Driven Control for Linear Discrete-Time Systems with Time-Varying Delays
%%author:Pengfei Kong from ECUST

% Discrete-time state-space model matrices
A = [0.9 0.5; 0.8 0.1];
B = [1 ;0.5];
d_max = 5;
d_min = 1;
n = size(A, 1);
m = size(B,2);
L = size(A,1) + 2*size(B,1) + d_max + 4;%data length

% State-feedback control gain
%K = [102.9100 80.7916];
Ad = [0.3 0;0.8 0.5];
% Simulation parameters
x0 = [2; -2];  
N = 25;  % Number of time steps

% Initialize state and input histories
x_hist = zeros(2, N);
u_hist = zeros(1, N);
d_hist = zeros(1, N); % Array to store delay at each time step

% Set initial state
x_hist(:, 1) = x0;
u_hist = 2*rand(1, N) -1;
% Simulation loop
for k = 1:N-1
    % Generate a random delay between 1 and 5
    dk = randi([d_min, d_max]);
    d_hist(k) = dk; % Record the delay
    
    if k > dk
        % State feedback control with random delay
         x_hist(:, k + 1) = A * x_hist(:, k) + B * u_hist(k) + Ad * x_hist(:, k - dk);
    else
        % No control applied for the first dk steps
        x_hist(:, k + 1) = A * x_hist(:, k) + B * u_hist(k);
    end
end

% Plot results
figure(1);
plot(0:N-1, x_hist(1, :), 'b-', 'LineWidth', 1.2,'DisplayName', 'State 1');
hold on;
plot(0:N-1, x_hist(2, :), 'r-', 'LineWidth', 1.2,'DisplayName', 'State 2');
xlabel('Time step');
ylabel('State Value');
title('State Trajectories');
legend;
hold off;


%data matrix
X_0 = zeros(n,L);
X_1 = zeros(n,L);
X_d = zeros(n,L);
U_0 = zeros(m,L);

X_0 = x_hist(:,6:6+L-1);
X_1 = x_hist(:,7:6+L);
U_0 = u_hist(:,6:6+L-1);
for k=1:15
    X_d(:,k) = x_hist(:,k+5-d_hist(k+5));
end

%LMI
S1 = sdpvar(n, n, 'symmetric');
S2 = sdpvar(n, n, 'full');
S3 = sdpvar(n, n, 'full');
Q = sdpvar(n, n, 'symmetric');
Y = sdpvar(n, n, 'symmetric');
Q1 = sdpvar(L, n, 'full');
Q2 = sdpvar(L, n, 'full');



Phi = [
    -S1, -(X_1*Q1)' + S2',  S2', zeros(n, n), S1', S1';
    -(X_1*Q1) + S2, S3 + S3', S3', -X_1*Q2, zeros(n, n), zeros(n, n);
    S2, S3, -S1, zeros(n, n), zeros(n, n), zeros(n, n);
    zeros(n, n), -(X_1*Q2)', zeros(n, n), -Q, zeros(n, n), zeros(n, n);
    S1, zeros(n, n), zeros(n, n), zeros(n, n), -Q, zeros(n, n);
    S1, zeros(n, n),zeros(n, n), zeros(n, n), zeros(n, n), -(1/(d_max - d_min))*Y
    ];

% 添加LMI约束
F = [Phi <= 0,  S1 >= 0, Q >= 0, Y >= 0 , Y <= Q ,...
     X_d *Q1 == X_0*Q2,X_0*Q2 == 0 ,S1 == X_0*Q1 , Q == X_d*Q2];


options = sdpsettings('solver', 'mosek', 'verbose', 3);
diagnostics = optimize(F, [], options)


if diagnostics.problem == 0
    disp('LMI problem is feasible.');
    Q1_value = value(Q1);

    K = U_0 * Q1_value /(value(S1)); 
    disp('The control gain K is:');
    disp(K);
else
    error('LMI problem is not feasible.');
end



for k = 1:N-1
    % Generate a random delay between 1 and 5
    dk = randi([1, 5]);
    d_hist2(k) = dk; % Record the delay
    u_hist(k) = K * x_hist(:, k);
    if k > dk
        % State feedback control with random delay
         x_hist(:, k + 1) = A * x_hist(:, k) + B * u_hist(k) + Ad * x_hist(:, k - dk);
    else
        % No control applied for the first dk steps
        x_hist(:, k + 1) = A * x_hist(:, k) + B * u_hist(k);
    end
end

% Plot results
figure(2);
plot(0:N-1, x_hist(1, :), 'b-', 'LineWidth', 1.2,'DisplayName', 'State 1');
hold on;
plot(0:N-1, x_hist(2, :), 'r-','LineWidth', 1.2, 'DisplayName', 'State 2');
xlabel('Time step');
ylabel('State Value');
title('State Trajectories');
xlim([0, N-1]);
legend;

hold off;


%delay figure
figure;
stem(0:1:24, d_hist(1:25), 'filled', 'LineWidth', 1.2);

xlabel('Time step', 'FontSize', 12);    
ylabel('Delay value', 'FontSize', 12);  