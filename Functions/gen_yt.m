

M = 4;
N = 8;
shape = "rrc";
alpha = 0.4;
Q = 8;
T = 1 / 15000;
Ts = T / M;
res = 10;
dt = Ts / res;
t_range = 0:dt:(N*T + Q*Ts);
N0 = 0;
z_t = sqrt(N0/2) * (randn(1,length(t_range)) + 1j*randn(1,length(t_range)));
Phi_i = 1;
tau_i = 0;
v_i = 0;
xDD = [1;zeros(N*M-1,1)];

% function [y_t,yDD] = gen_yt(z_t,Phi_i,tau_i,v_i,xDD,N,M,T,shape,alpha,Q,res)

% Parameters
Ts = T / M;
F0 = 1 / (N*T);
dt = Ts / res;
t_range = 0:dt:(N*T + Q*Ts);

% Generate time-domain TX signal
[x_t,g_t] = gen_xt(xDD,N,M,T,shape,alpha,Q,res);

% Generate channel signal
r_t_noiseless = zeros(1,length(t_range));
for i = 1:length(Phi_i)

    t_new = t_range - tau_i(i);
    x_t_shift = interp1(t_range, x_t, t_new, 'spline', 'extrap');

    exp_t = exp(1j.*2.*pi.*v_i(i).*t_new);

    r_t_noiseless = r_t_noiseless + Phi_i(i) .* x_t_shift .* exp_t;

end

% Add noise to RX signal
r_t = r_t_noiseless + z_t;

% Get normal value for pulses
a_t_sample = gen_pulse(0:dt:Q*Ts,shape,Ts,Q,alpha);
norm_val = sqrt(sum(abs(a_t_sample).^2) * dt) * sqrt(N);

y_t = zeros(length(t_range),N);
yDD = zeros(M,N);
for m = 0:M-1
    for n = 0:N-1

        % Set up delay and Doppler taps
        t = m*Ts;
        f = n*F0;

        % Set up exponential function
        exp_t = exp(-1j.*2.*pi.*f.*(t_range-t));

        % Set up TX pulse train p*(-t)
        n2 = (-1:N).';
        t_range_at = (t_range-t) - n2 * T;
        a_t = gen_pulse(t_range_at,shape,Ts,Q,alpha);
        p_t = sum(a_t, 1) / norm_val;

        y_int_fun = r_t .* conj(p_t) .* exp_t;

        yDD(m+1,n+1) = sum(abs(y_int_fun).^2) * dt;

        % figure(1)
        % plot(t_range_postconv/Ts,abs(y_t(:,n+1)));
        % xlim([0 N*T/Ts])
        % 1;


        % t_taps = (0:(M-1))*Ts;
        % t_differences = t_range - t_taps.';
        % loc = zeros(M,1);
        % for m = 0:M-1
        %     [~,loc(m+1)] = min(abs(t_differences(m+1,:)));
        % end
        % 
        % yDD(:,n+1) = y_t(loc,n+1).';
        % 
        % 1;


    end
end

yDD = yDD(:);

% end
