% Square-wave response with DC offsets
% Shows steady-state convergence levels

clear; clc; close all;

A = 1e-10; B = 1e-4;
sigma_on = 0.45; sigma_off = 0.013;
x_on = 0.06; x_off = 0.4;
sigma_p = 4e-5; beta = 500;
Gm = 0.025; b = 4.7; a = 7.2e-6;

V = 0.45;
f = 5.8e6;
T = 1/f;

Vo_list = [-0.0428 0.0086 0.0201 0.0400 0.0662 0.0957 0.1267 0.1581 0.1892];

G = @(x,v) Gm*x + a*exp(b*sqrt(abs(v))).*(1-x);
g_plus  = @(x,v) B*sinh(v./sigma_on).*exp(-(x.^2)/(x_on^2)).*exp((G(x,v).*v.^2)/sigma_p);
g_minus = @(x,v) A*sinh(v./sigma_off).*exp(-(x_off^2)./(x.^2)).*exp(-(1 + beta*G(x,v).*v.^2));

figure; hold on;

for Vo = Vo_list

    x = 1;
    t_hist = 1e-10;
    x_hist = x;

    for n = 1:round(100/T)

        x = x + g_plus(x,  V+Vo)*(T/2);
        x = x + g_minus(x, -V+Vo)*(T/2);
        x = min(max(x,0),1);

        if mod(n,50)==0
            t_hist(end+1) = n*T;
            x_hist(end+1) = x;
        end
    end

    semilogx(t_hist, x_hist, 'LineWidth', 2);
end

xlabel('t / s'); ylabel('x');
title('Memristor state response (square wave)');
grid on; xlim([1e-10 100]); ylim([0 1.05]);
legend(string(Vo_list),'Location','eastoutside');
