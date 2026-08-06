// =====================================================================
// Q7. Hump shape in capital  --  derived from soe1.mod
//     Model, steady state and shocks blocks are byte-identical to
//     soe1.mod. Only the post-estimation section is new.
// =====================================================================

// Standard Open Economy DSGE Model
// Based on Uribe and Schmitt-Grohé chapter 4
// by Carlos Goncalves

var c, w, lambda, k, d, h, r, A, y, ca, i, tb, tb_y, ca_y;
varexo eps_A;

parameters beta, alpha, delta, omega, sigma, phi, rstar, psi, rho_A, sigma_A, d_bar;

// Parameter values
beta  = 0.96;    // Discount factor
alpha = 0.32;    // Capital share
delta = 0.1;     // Depreciation rate
omega = 1.5;    // Labor disutility parameter  
sigma = 2;       // Risk aversion
phi   = 0.01;     // Capital adjustment cost parameter
rstar = (1/beta)-1;    // International interest rate
psi   = 0.01;   // Debt elasticity parameter
rho_A = 0.4;    // Productivity persistence
sigma_A = 0.01;  // Productivity shock std
d_bar = 0.7;     // Steady state debt level

model;

// MgU of consumption = shadow price 
lambda=c-(h^omega)/omega;

// Choosing debt gives Euler
lambda = beta*(1+r(+1))*lambda(+1);

// But you can also save in capital, choosing K(t+1) 
(1+phi*(k(0)-k(-1)))*lambda = beta*lambda(+1)*(alpha*A(+1)*k(0)^(alpha-1)*h(+1)^(1-alpha) + (1-delta)+ phi*(k(+1)-k(0)));

// Investment including 
i = k - (1-delta)*k(-1);

// MgDisutility of labor = MgU of consumption*real wages 

w = h^(omega-1);

// demand for labor from firms

w= (1-alpha)*A*(k(-1)^alpha)*h^(-alpha);

// Production function

y = A*(k(-1)^alpha)*(h^(1-alpha));

// Budget constraint : to increase capital by i, you have to spend i+adj.costs

d = c + (i + 1/2*phi*(k(0)-k(-1))^2) + (1+r)*d(-1) - y;

// interest rate (debt-elastic around d_bar)
r(+1) = rstar + psi*(exp(d-d_bar) - 1);

// current account and trade balance-- definition

ca=-d+d(-1);

tb-r*d = ca;

tb_y=tb/y;

ca_y = ca/y;

// Productivity process
log(A) = rho_A*log(A(-1)) + eps_A;

end;

initval;
A = 1;
h = 1.1;
k = 4;
y = 1.7;
i = 0.4;
d = d_bar;
r = rstar;
c = 1.2;
ca=0;
tb=0.02;
end;

steady_state_model;
h = ((1-alpha)*((1/beta - (1-delta)) / alpha)^(alpha/(alpha-1)))^(1/(omega -1)); 
k = h*((1/beta - (1-delta)) / alpha)^(1/(alpha-1)) ;
y =(k^alpha)*(h^(1-alpha));
c= y-delta*k-rstar*d_bar;
i = delta*k;
w= (1-alpha)*h^(-alpha)*k^alpha;
lambda = c-(h^omega)/omega;
r = rstar ;
d=d_bar;
ca=0;
A=1;
tb = rstar*d_bar; 
tb_y=tb/y;
ca_y=0;
end;

steady;  

shocks;
var eps_A; stderr sigma_A;
end;


// =====================================================================
// Q7. Hump shape in capital
//
// Everything above this line is soe1.mod, unchanged.
// The only thing that changes is the capital adjustment cost phi, which
// is set inside the loop with set_param_value - the parameter block above
// is untouched.
//
// Why it works: the adjustment cost is (phi/2)*(k_t - k_{t-1})^2, so the
// marginal cost of moving capital in one go rises with phi. With phi
// almost zero the firm jumps straight to the new desired capital stock in
// the period of the shock, and k decays monotonically afterwards. Make phi
// big and it is cheaper to spread the same investment over several
// periods: k keeps rising for a while before TFP has decayed enough to
// turn it around. That interior peak is the hump.
//
// order=1 rather than the original order=2: at second order Dynare gets
// the IRFs by Monte Carlo and they come out jagged, which hides the shape.
// irf=40 rather than 10 because the hump peaks well after period 10.
// =====================================================================


% Figures go to output/figures/, next to the ones from R. The .mod files are
% run from inside dynare/, so this is one level up.
outdir = fullfile('..','output','figures');

phi_values = [0.01 0.5 2 6];
colors = {'blue', 'red', 'green', 'magenta'};
markers = {'o', 's', '^', 'd'};

figure('Name','Q7: capital IRF and the adjustment cost','Position',[100 100 1100 700]);
set(gcf,'Color','white'); hold on;

k_pos = strmatch('k', M_.endo_names, 'exact');

for j = 1:length(phi_values)
    set_param_value('phi', phi_values(j));
    steady;
    stoch_simul(order=1, irf=40, periods=0, nograph);

    k_irf = 100 * oo_.irfs.k_eps_A / oo_.steady_state(k_pos);
    periods = 0:options_.irf-1;

    plot(periods, k_irf, 'Color', colors{j}, 'Marker', markers{j}, ...
         'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', colors{j}, ...
         'DisplayName', ['\phi = ' num2str(phi_values(j))]);

    [peak, ipeak] = max(k_irf);
    fprintf('phi = %5.2f : peak %+7.4f%% at period %2d\n', ...
            phi_values(j), peak, ipeak-1);
end

plot([0 options_.irf-1],[0 0],'--','LineWidth',1.5,'Color',[0.5 0.5 0.5], ...
     'HandleVisibility','off');
xlabel('Periods after shock','FontSize',13,'FontWeight','bold');
ylabel('% deviation from steady state','FontSize',13,'FontWeight','bold');
title('Capital response to a TFP shock: raising \phi produces the hump', ...
      'FontSize',15,'FontWeight','bold');
legend('Location','best','FontSize',12); grid on; grid minor;
xlim([0 options_.irf-1]); hold off;
print('-dpng','-r300',fullfile(outdir,'q7_capital_hump.png'));

// The question says "adjustment cost OR rho". Same experiment with rho_A,
// holding phi at its original value.
set_param_value('phi', 0.01);
rho_values = [0.4 0.7 0.9 0.95];

figure('Name','Q7: capital IRF and TFP persistence','Position',[150 150 1100 700]);
set(gcf,'Color','white'); hold on;

for j = 1:length(rho_values)
    set_param_value('rho_A', rho_values(j));
    steady;
    stoch_simul(order=1, irf=40, periods=0, nograph);

    k_irf = 100 * oo_.irfs.k_eps_A / oo_.steady_state(k_pos);
    plot(0:options_.irf-1, k_irf, 'Color', colors{j}, 'Marker', markers{j}, ...
         'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerFaceColor', colors{j}, ...
         'DisplayName', ['\rho_A = ' num2str(rho_values(j))]);

    [peak, ipeak] = max(k_irf);
    fprintf('rho_A = %4.2f : peak %+7.4f%% at period %2d\n', ...
            rho_values(j), peak, ipeak-1);
end

plot([0 options_.irf-1],[0 0],'--','LineWidth',1.5,'Color',[0.5 0.5 0.5], ...
     'HandleVisibility','off');
xlabel('Periods after shock','FontSize',13,'FontWeight','bold');
ylabel('% deviation from steady state','FontSize',13,'FontWeight','bold');
title('Capital response to a TFP shock: raising \rho_A also delays the peak', ...
      'FontSize',15,'FontWeight','bold');
legend('Location','best','FontSize',12); grid on; grid minor;
xlim([0 options_.irf-1]); hold off;
print('-dpng','-r300',fullfile(outdir,'q7_capital_hump_rho.png'));

fprintf('\nA peak at period 0 is no hump. A peak at period 2 or later is.\n');
