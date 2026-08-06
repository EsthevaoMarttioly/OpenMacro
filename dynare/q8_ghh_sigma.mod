// =====================================================================
// Q8. sigma and the GHH aggregator  --  derived from
//     GPU_2010_byCESG_v110.mod (Garcia-Cicco, Pancrazi, Uribe 2010).
//     Three lines added, all marked "Q8"; everything else is the
//     original file down to write_latex_dynamic_model.
// =====================================================================

/*
 * This file replicates the model studied in:
 * García-Cicco, Javier and Pancrazi, Roberto and Uribe, Martín (2010): "Real Business Cycles
 * in Emerging Countries", American Economic Review, 100(5), pp. 2510-2531.
 */

// Dynare options to avoid filesystem conflicts
@#define USE_DLL = 0
@#define bytecode = 0

@#define RBC =0 
//set to 1 for RBC model and to 0 for Financial Frictions Model

var c       $c$ 
    k       $k$ 
    a       $a$ 
    h       $h$ 
    d       $d$ 
    y       $y$ 
    invest  $i$  
    tb      $tb$ 
    mu_c    ${MU_C}$ 
    tb_y    ${\frac{TB}{Y}}$ 
    growth_y     ${\Delta Y}$
    growth_c     ${\Delta C}$
    growth_invest ${\Delta I}$
    g       ${g}$
    r       ${r}$
    % financial shock added to compete with Gita's paper
    mu      ${\mu}$
    % preference shock added to compete with Gita's paper; 
    % whenever I see this shock I think my dear Bernardo Guimaraes would be in shock (pun intended) 
    nu      ${\nu}$
    G       ${G}$          % Q8: GHH consumption aggregator
    @#if RBC == 0
    % this is a government demand shock    
    s       ${s}$
    @# endif
; 

predetermined_variables k d;

%Define parameters
parameters beta     ${\beta}$ 
        gamma       ${\gamma}$ 
        delta       ${\delta}$
        alpha       ${\alpha}$
        psi         ${\psi}$
        omega       ${\omega}$
        theta       ${\theta}$
        phi         ${\phi}$
        dbar        ${\bar d}$
        gbar        ${\bar g}$
        rho_a       ${\rho_a}$
        rho_g       ${\rho_g}$
        rho_nu      ${\rho_\nu}$ % governing U(.) shock
        rho_mu      ${\rho_\mu}$ % governing financial shock
        rho_s       ${\rho_s}$ % gov demand push
    @#if RBC == 0
        s_share     ${sshare}$
        S           ${S}$
    @# endif
;

varexo eps_a ${\varepsilon_a}$
        eps_g ${\varepsilon_g}$ 
        eps_nu ${\varepsilon_\nu}$
        eps_mu ${\varepsilon_\mu}$
    @#if RBC == 0
        eps_s ${\varepsilon_s}$
    @# endif
;

% this is the calibration if we are NOT adding stuff to Gita's paper        
@#if RBC == 1
    gbar  = 1.0050; %Gross long term growth rate
    rho_g = 0.8280; %Serial correlation of innovation in permanent technology shock
    rho_a = 0.7650; %Serial correlation of transitory technology shock
    phi   = 3.3000; %Adjustment cost parameter
@# else
    gbar  = 1.009890776104921; % this is weird, why gbar would be different?
    rho_g = 0.323027844166870; % when AG and Uribe et al estimate these models, they get very different results for this
    rho_a = 0.864571930755821;
    phi   = 4.810804146604144;
@# endif
            
rho_nu = 0.850328786147732;
rho_s  = 0.205034667802314;
rho_mu = 0.906802888826967;

%From Table 2
gamma = 2; %intertemporal elasticity of substitution
delta = 1.03^4-1;%0.03; %Depreciation rate
alpha = 0.32; %Capital elasticity of the production function
omega = 1.6; %exponent of labor in utility function
theta = 1.4*omega;
beta = 0.98^4; % discount factor, calibrating to annual -- this number is a bit off
dbar = 0.007;

@#if RBC == 1
    psi = 0.001; % seems low for an EM
@# else
%psi is sooo different: quizz to you: why???    
psi = 2.867166241970346; %parameter governing the debt elasticity of the interest rate.
s_share = 0.10; %Share of public spending in GDP
@# endif
        
        
model;
#RSTAR = 1/beta * gbar^gamma; %World interest rate
% this # means it is a pre-computed variable, before the model is run, we can use it throughout

%1. Interest Rate
r = RSTAR + psi*(exp(d-dbar) - 1)+exp(mu-1)-1;

%2. Marginal utility of consumption
mu_c = nu * (c - theta/omega*h^omega)^(-gamma);

%2b. Q8: the GHH aggregator itself, so we can plot it
G = c - theta/omega*h^omega;

%3. Resource constraint (see the remark on the fixed typo in the preamble)
@#if RBC == 1
    y= log(tb) + c + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# else
% remember, here we have government !
    y= log(tb) + c + s + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# endif

%4. Trade balance
% from this equation, you can infer that g is a gross rate ...
log(tb)= d - d(+1)*g/r;

%5. Definition output
y= a*k^alpha*(g*h)^(1-alpha);

%6. Definition investment
invest= k(+1)*g - (1-delta) *k;

%7. Euler equation
% from here you should be able to realize that r is a GROSS rate
mu_c= beta/g^gamma*r*mu_c(+1);

%8. First order condition labor
theta*h^(omega-1)=(1-alpha)*a*g^(1-alpha)*(k/h)^alpha;

%9. First order condition investment
mu_c*(1+phi*(k(+1)/k*g-gbar))= beta/g^gamma*mu_c(+1)*(1-delta+alpha*a(+1)*(g(+1)*h(+1)/k(+1))^(1-alpha) +phi*k(+2)/k(+1)*g(+1)*(k(+2)/k(+1)*g(+1)-gbar) - phi/2*(k(+2)/k(+1)*g(+1)-gbar)^2);

%10. Definition trade-balance to output ratio
log(tb_y) = log(tb)/y; 

%11. Output growth
% pay attention here: y is deviation from trend only, so you have to multiply by g
growth_y= y/y(-1)*g(-1);

%12. Consumption growth
growth_c = c/c(-1)*g(-1);

%13. Investment growth
growth_invest = invest/invest(-1)*g(-1);

%14. Law of Motion = temporary TFP
log(a)=rho_a * log(a(-1))+eps_a; 

%15. Law of Motion = TFP Growth
log(g/gbar)=rho_g*log(g(-1)/gbar)+eps_g; 

%16. Law of Motion = Preference shock
log(nu) =rho_nu * log(nu(-1))+eps_nu;

%17. Law of Motion = exogenous stochastic country premium shock
log(mu)= rho_mu * log(mu(-1))+eps_mu;

@#if RBC == 1
    
@# else
    %18. Law of Motion = Exogenous spending shock
    log(s/S)= rho_s * log(s(-1)/S) + eps_s;
@# endif

end;

steady_state_model;
    r   = 1/beta*gbar^gamma; % World interest rate
    d   = dbar; %foreign debt
    % pen and paper get you here:
    k_over_gh =((gbar^gamma/beta-1+delta)/alpha)^(1/(alpha-1)); %k/(g*h)
    % pen and paper get you here:
    h   = ((1-alpha)*gbar*k_over_gh^alpha/theta)^(1/(omega-1)); %hours
    % definition get you here:
    k   = k_over_gh*gbar*h; %capital
    invest = (gbar-1+delta)*k; %investment
    y   = k^alpha*(h*gbar)^(1-alpha); %output
    @#if RBC == 1
        s = 0;    
    @# else
        s = y*s_share;
        S = s;
    @# endif
    c   = (gbar/r-1)*d +y-s-invest; %Consumption
    tb  = y - c - s - invest; %Trade balance
    tb_y = tb /y;
    mu_c = (c - theta/omega*h^omega)^(-gamma); %marginal utility of wealth
    G    = c - theta/omega*h^omega; % Q8
    a   = 1; %productivity shock 
    g   = gbar; %Growth rate of nonstationary productivity shock
    growth_c = g;
    growth_invest = g;
    growth_y = g;
    nu  = 1; % why not 0 ? They are multiplicative, go back to the equations !
    mu  = 1;
    tb  = exp(tb); 
    % explain this on black board: when tb=1, it is really (trade balance)/(trade balance (SS)) = 1
    % log(1) = 0 
    tb_y= exp(tb_y);
end;

shocks;
@#if RBC == 1
    var eps_a; stderr 0.0270;
    var eps_g; stderr 0.0300;
    var eps_nu; stderr  0;
    var eps_mu; stderr  0;
@# else
    var eps_a; stderr 0.033055089525252;
    var eps_g; stderr 0.010561526060797;
    var eps_nu; stderr  0.539099453618175;
    var eps_s; stderr   0.018834174505537;
    var eps_mu; stderr  0.057195449717680;
@# endif
end;


// =====================================================================
// Q8. What sigma (here: gamma) does to the GHH aggregator G
//
// Everything above is GPU_2010_byCESG_v110.mod with exactly three lines
// added, each marked "Q8": G in the var block, its definition in the
// model block, and its steady state. gamma itself is never edited in the
// parameter block - the loop below sets it with set_param_value.
//
// Where sigma lives in this model. Preferences are GHH:
//
//     U = G^(1-gamma)/(1-gamma),      G = c - (theta/omega) h^omega
//
// so the file's "gamma" is the sigma of the question, and G is the whole
// argument of the utility function. Equation 2 of the model,
//
//     mu_c = nu * (c - theta/omega*h^omega)^(-gamma)
//
// is the marginal utility of G. The point of GHH is that the labour FOC,
// theta*h^(omega-1) = MPL, has no consumption in it: hours depend only on
// the wage, and the household smooths G, not c. So plotting c alone
// hides the object the model is actually about, which is why the
// question asks for the whole G.
//
// What to expect. gamma is the inverse intertemporal elasticity, and
// 1/gamma is how willing the household is to move G across time. Push
// gamma to 5 or 6 and it is very unwilling: G becomes flat and the
// adjustment is thrown onto the trade balance instead. Take gamma to 1,
// which is log utility of G, and G moves a lot more.
// =====================================================================


% Figures go to output/figures/, next to the ones from R. The .mod files are
% run from inside dynare/, so this is one level up.
outdir = fullfile('..','output','figures');

gamma_values = [1 2 5 6];
colors  = {'blue', 'red', 'green', 'magenta'};
markers = {'o', 's', '^', 'd'};

shock_list = {'eps_a', 'eps_g'};
shock_name = {'transitory TFP (\epsilon_a)', 'trend TFP (\epsilon_g)'};

sigma_ratio = zeros(length(gamma_values), 2);

for s = 1:2
    % Keep explicit handles. stoch_simul runs between the plot calls and can
    % leave gcf pointing somewhere else, which makes print() fail with
    % "the figure is invalid".
    fh = figure('Name',['Q8: G and c, ' shock_list{s}],'Position',[100+40*s 100 1200 520]);
    set(fh,'Color','white');
    ax1 = subplot(1,2,1,'Parent',fh); hold(ax1,'on');
    ax2 = subplot(1,2,2,'Parent',fh); hold(ax2,'on');

    for j = 1:length(gamma_values)
        set_param_value('gamma', gamma_values(j));
        steady;
        stoch_simul(loglinear, order=1, irf=20, nograph) G c y tb_y;

        % loglinear => IRFs are already percent deviations
        G_irf = eval(['oo_.irfs.G_' shock_list{s}]);
        c_irf = eval(['oo_.irfs.c_' shock_list{s}]);

        plot(ax1, 0:options_.irf-1, 100*G_irf, 'Color', colors{j}, 'Marker', markers{j}, ...
             'LineWidth', 2.4, 'MarkerSize', 6, 'MarkerFaceColor', colors{j}, ...
             'DisplayName', ['\gamma = ' num2str(gamma_values(j))]);

        plot(ax2, 0:options_.irf-1, 100*c_irf, 'Color', colors{j}, 'Marker', markers{j}, ...
             'LineWidth', 2.4, 'MarkerSize', 6, 'MarkerFaceColor', colors{j}, ...
             'DisplayName', ['\gamma = ' num2str(gamma_values(j))]);

        if s == 1
            iy = strmatch('y', oo_.var_list, 'exact');
            ic = strmatch('c', oo_.var_list, 'exact');
            iG = strmatch('G', oo_.var_list, 'exact');
            sigma_ratio(j,1) = sqrt(oo_.var(ic,ic)/oo_.var(iy,iy));
            sigma_ratio(j,2) = sqrt(oo_.var(iG,iG)/oo_.var(iy,iy));
        end
    end

    for ax = [ax1 ax2]
        plot(ax, [0 options_.irf-1],[0 0],'--','LineWidth',1.2,'Color',[0.5 0.5 0.5], ...
             'HandleVisibility','off');
        xlabel(ax,'Periods after shock','FontSize',12,'FontWeight','bold');
        ylabel(ax,'% deviation from SS','FontSize',12,'FontWeight','bold');
        grid(ax,'on'); grid(ax,'minor'); xlim(ax,[0 options_.irf-1]);
        legend(ax,'Location','best','FontSize',10);
        hold(ax,'off');
    end
    title(ax1, ['G = c - (\theta/\omega)h^\omega  --  ' shock_name{s}], ...
          'FontSize',13,'FontWeight','bold');
    title(ax2, ['Consumption only  --  ' shock_name{s}], ...
          'FontSize',13,'FontWeight','bold');

    print(fh,'-dpng','-r300',fullfile(outdir,['q8_G_vs_c_' shock_list{s} '.png']));
end

fprintf('\n=== Q8: volatility ratios, loglinear, all shocks on ===\n');
fprintf('gamma   sigma_C/sigma_Y   sigma_G/sigma_Y\n');
fprintf('-----   ---------------   ---------------\n');
for j = 1:length(gamma_values)
    fprintf('%5.2f   %15.3f   %15.3f\n', gamma_values(j), sigma_ratio(j,1), sigma_ratio(j,2));
end
fprintf('\ngamma = 1 is log utility of G. Compare the two columns: c and G\n');
fprintf('are not the same object, and the gap is the labour term.\n');
