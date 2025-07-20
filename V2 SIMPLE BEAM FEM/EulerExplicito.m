function [x, v, a, t] = EulerExplicito(x_0, v_0, a_0, N_step, delta_t, M, K, mu)

x = zeros(N_step, 1);
v = zeros(N_step, 1);
a = zeros(N_step, 1);
t = zeros(N_step, 1);

a(1) = -K/M*x_0;
v(1) = v_0 + a(1)*delta_t;
x(1) = x_0 + v(1)*delta_t;


for i = 2:N_step
    a(i) = -K/M*x(i - 1) - mu*v(i - 1);
    v(i) = v(i - 1) + a(i)*delta_t;
    x(i) = x(i - 1) + v(i)*delta_t;
    t(i) = t(i- 1) + delta_t;
end

end