// =====================================================================
// Q9. rho close to one  --  derived from soe1.mod
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
// Q9. rho close to one, and what happens to sigma_C / sigma_Y
//
// Everything above this line is soe1.mod, unchanged. Only rho_A moves,
// through set_param_value inside the loop.
//
// What to look for:
//   1. sigma_C/sigma_Y rises as rho_A goes to 1. A near-permanent TFP
//      shock raises permanent income roughly one for one with current
//      output, so consumption moves at least as much as output. That is
//      the AG mechanism: excess consumption volatility does not need
//      credit constraints, it needs persistent shocks.
//   2. The trade balance turns countercyclical. With rho_A near 1 the
//      country borrows to invest and to consume, so tb/y falls on impact.
//   3. The tb/y path at rho_A = 0.95 looks like the one at rho_A = 0.05
//      on impact, but it comes back to zero far more slowly. That is the
//      "volta bem devagar" in the question - near a unit root the only
//      thing pulling debt back is the debt-elastic premium psi, which is
//      tiny (0.01), so the half-life is enormous.
//
// hp_filter=100 so the moments are comparable to the annual HP-filtered
// moments in the empirical part; sigma ratios are taken as coefficients
// of variation, which to first order equal the ratio of log standard
// deviations. irf=60 to see the slow return.
// =====================================================================


% Figures go to output/figures/, next to the ones from R. The .mod files are
% run from inside dynare/, so this is one level up.
outdir = fullfile('..','output','figures');

rho_values = [0.05 0.4 0.7 0.9 0.95];
colors  = {'blue', 'red', 'green', 'magenta', 'black'};
markers = {'o', 's', '^', 'd', 'p'};

results = zeros(length(rho_values), 5);
tb_paths = zeros(length(rho_values), 60);

y_pos = strmatch('y',    M_.endo_names, 'exact');
c_pos = strmatch('c',    M_.endo_names, 'exact');
t_pos = strmatch('tb_y', M_.endo_names, 'exact');

for j = 1:length(rho_values)
    set_param_value('rho_A', rho_values(j));
    steady;
    stoch_simul(order=1, hp_filter=100, irf=60, periods=0, nograph) y c tb_y;

    iy = strmatch('y',    oo_.var_list, 'exact');
    ic = strmatch('c',    oo_.var_list, 'exact');
    it = strmatch('tb_y', oo_.var_list, 'exact');

    // coefficients of variation: sd relative to the steady-state level
    sd_y = sqrt(oo_.var(iy,iy)) / oo_.steady_state(y_pos);
    sd_c = sqrt(oo_.var(ic,ic)) / oo_.steady_state(c_pos);

    // correlation of the trade balance ratio with output
    corr_tb_y = oo_.var(it,iy) / sqrt(oo_.var(it,it)*oo_.var(iy,iy));

    results(j,:) = [rho_values(j), 100*sd_y, 100*sd_c, sd_c/sd_y, corr_tb_y];
    tb_paths(j,:) = oo_.irfs.tb_y_eps_A;
end

fprintf('\n=== Q9: sigma_C / sigma_Y as rho_A goes to one ===\n');
fprintf('rho_A   sigma_y(%%)  sigma_c(%%)  sigma_c/sigma_y   corr(tb/y, y)\n');
fprintf('-----   ----------  ----------  ---------------   -------------\n');
for j = 1:length(rho_values)
    fprintf('%5.2f   %10.3f  %10.3f  %15.3f   %13.3f\n', results(j,:));
end

% Explicit handles: gcf can be left pointing elsewhere after stoch_simul, and
% print() then fails with "the figure is invalid".
fh1 = figure('Name','Q9: trade balance response, low vs high persistence','Position',[100 100 1100 700]);
set(fh1,'Color','white'); hold on;
for j = 1:length(rho_values)
    plot(0:59, tb_paths(j,:), 'Color', colors{j}, 'Marker', markers{j}, ...
         'LineWidth', 2.2, 'MarkerSize', 5, 'MarkerFaceColor', colors{j}, ...
         'MarkerIndices', 1:4:60, ...
         'DisplayName', ['\rho_A = ' num2str(rho_values(j))]);
end
plot([0 59],[0 0],'--','LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
xlabel('Periods after shock','FontSize',13,'FontWeight','bold');
ylabel('tb/y, deviation from SS','FontSize',13,'FontWeight','bold');
title('Trade balance: rho_A near 1 starts like rho_A near 0 but returns very slowly', ...
      'FontSize',14,'FontWeight','bold');
legend('Location','best','FontSize',12); grid on; grid minor;
xlim([0 59]); hold off;
print(fh1,'-dpng','-r300',fullfile(outdir,'q9_tb_persistence.png'));

fh2 = figure('Name','Q9: sigma_C / sigma_Y','Position',[150 150 900 600]);
set(fh2,'Color','white');
plot(results(:,1), results(:,4), '-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'Color', 'blue', 'MarkerFaceColor', 'blue');
hold on;
plot([min(rho_values) max(rho_values)], [1 1], '--', 'LineWidth', 1.5, ...
     'Color', [0.5 0.5 0.5]);
xlabel('\rho_A','FontSize',13,'FontWeight','bold');
ylabel('\sigma_C / \sigma_Y','FontSize',13,'FontWeight','bold');
title('Relative consumption volatility against TFP persistence', ...
      'FontSize',14,'FontWeight','bold');
grid on; grid minor; hold off;
print(fh2,'-dpng','-r300',fullfile(outdir,'q9_sigma_ratio.png'));
