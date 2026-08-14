// =========================================================================
// DSGE_uribe_yue_IRFmatch_Brazil.mod
//
// ****************************************************************
// *** BRAZIL ONLY. Sample: 1996:Q1 - 2026:Q1 (single-country VAR ***
// *** replication, uribe_yue_VAR_Brazil.mlx / Table 6.1.BR).     ***
// *** NOT the paper's original 7-country panel (1994-2001) --     ***
// *** every result below reflects Brazil's own, much longer and   ***
// *** more current, estimated dynamics.                           ***
// ****************************************************************
//
// Small open economy model with working-capital-in-advance, external habit
// formation, debt-elastic country risk premium, and time-to-build capital.
// Uribe and Yue, "Country Spreads and Emerging Countries: Who Drives Whom?"
// (= Chapter 6 of Schmitt-Grohe and Uribe, "Open Economy Macroeconomics").
//
// THIS VERSION ESTIMATES {PSSI, PHI, ETA, MU} BY IRF MATCHING -- the same
// four "deep" adjustment-cost / friction parameters the paper itself
// calibrates (Section 3.1.4, eq. 17) by matching the model's theoretical
// impulse responses against the estimated panel-VAR's impulse responses.
// Following the paper's methodology exactly: 4 variables (output,
// investment, trade balance/GDP, country rate -- NOT the world rate) x 2
// shocks (country-spread, US rate) x 24 quarters = 192 target points, each
// shock sized at 0.01 (matching Figures 2/3's "unit innovation" / "one
// percentage point increase"). The target IRFs are produced by
// uribe_yue_VAR_Brazil_target_export.m (run after Section 3 of
// uribe_yue_VAR_Brazil.mlx) and loaded here from
// IRF_target_uribe_yue_Brazil.mat via matched_irfs; +
// method_of_moments(mom_method=IRF_MATCHING, ...)
// Equal-weighted across all 192
// points -- see the note above the matched_irfs; block for how this
// differs from the paper's inverse-variance weighting.
//
// WHAT'S DIFFERENT FROM THE BASELINE (DSGE_uribe_yue_IRFmatch.mod):
// (1) the ten coefficients of the two exogenous interest-rate processes,
// eqs. (6.20) and (6.22) -- ARRUS, ARRUS1, ARY, ARY1, ARIV, ARIV1, ARTBY,
// ARTBY1, ARR1, ARUSRUS1 -- are replaced with Table 6.1.BR's Brazil-only
// estimates (see the "Estimated VAR coefficients" section below); (2) the
// {PSSI,PHI,ETA,MU} STARTING VALUES (estimated_params_init;) are no
// longer the paper's own point -- Brazil's more persistent R/Rus process
// makes the model Blanchard-Kahn unstable there, so a stable point found
// by a dedicated prior-space scan is used instead (see the comment above
// estimated_params_init; for the full explanation). Priors themselves,
// and every other parameter (preferences, technology, steady-state
// targets), are unchanged from the baseline.
//
// Equation numbers in comments refer to Uribe and Yue's own numbering.
// =========================================================================

var
  c       // consumption
  h       // hours worked
  k       // capital stock
  d       // household debt / net foreign liabilities (level, can be negative)
  s0 s1 s2 s3   // investment projects at each stage of gestation (eq. 6.5-6.6)
  yy      // output
  tby     // trade-balance-to-GDP ratio (level, eq. 6.21)
  r       // country interest rate (log)
  rus     // world interest rate (log)
  iv      // aggregate investment (eq. 6.5)
  tb      // trade balance (level, eq. 6.21)
  la      // marginal utility of wealth (lambda)
  qq      // Tobin's Q
  nu0 nu1 nu2   // shadow prices of investment projects at stages 0-2
  er      // country-spread shock process (eq. 6.20)
  erus    // world-rate shock process (eq. 6.22)
;

varexo
  eps_r    // innovation to the country interest rate, eq. (6.20)
  eps_rus  // innovation to the world interest rate, eq. (6.22)
;

// Observables used by method_of_moments(mom_method=IRF_MATCHING) below --
// declared here, before model;
// .
varobs yy iv tby r;

// c, h, s0, and k are chosen one period in advance / are genuine states from
// the original toolkit's predetermined vector x, but Dynare's automatic
// classifier only recognizes a variable as predetermined when it has an
// EXPLICIT self-lag (-1) somewhere; none of these four ever get one (their
// only laws of motion are written as X(+1)=f(X,...), never X=f(X(-1),...)),
// so they must be declared explicitly to match the original state/control
// partition. 
predetermined_variables c h s0 k;

parameters
  BETTA GAMA DELTA ALFA PSSI ETA OMEGA PHI MU
  RSS RUSSS PREMIUM STB DBAR YYSS IVSS TBYSS
  ARUSRUS1 ARRUS ARRUS1 ARY ARY1 ARIV ARIV1 ARTBY ARTBY1 ARR1
  RHOR RHORUS
;

// ------------------------------------------------------------------
// Calibration (values as used by the authors to produce the paper's
// impulse-response figures; see the model's own steady-state / driver
// scripts). Nothing here is estimated inside Dynare -- these are the
// authors' point estimates / calibration targets, entered as constants.
// ------------------------------------------------------------------

// Preferences / technology
GAMA    = 2;             // risk-aversion (inverse IES)
OMEGA   = 1.455;         // labor-supply curvature (Frisch elast. = 1/(OMEGA-1))
DELTA   = 0.1/4;         // quarterly depreciation rate (10% annual)
ALFA    = 0.32;          // capital share

// Frictions (authors' calibrated/estimated values)
PSSI    = 0.000425;      // debt-adjustment cost intensity, Psi(d)=PSSI/2*(d-DBAR)^2
PHI     = 72.8268;       // capital-adjustment cost intensity, Phi(x)=x-PHI/2*(x-DELTA)^2
ETA     = 1.2023;        // fraction of wage bill subject to working-capital constraint
MU      = 0.2037;        // external habit intensity

// Interest rates / trade balance target
RUSSS   = 1+0.04/4;      // steady-state world interest rate (4% annual)
PREMIUM = 0.07/4;        // steady-state country spread (7% annual)
STB     = 0.02;          // steady-state trade-balance-to-GDP target
RSS     = RUSSS*(1+PREMIUM);
BETTA   = 1/RSS;         // discount factor consistent with steady-state bond Euler eq. (6.13)

// Estimated VAR coefficients for the country interest rate, eq. (6.20).
// BRAZIL ONLY (Table 6.1.BR, single-country 2SLS, N=119-120 for eqs.
// 1/2/3/5 and N=128 for the Rus AR(1); sample 1996:Q1-2026:Q1), NOT the
// original 7-country panel. Read directly off the printed Table 6.1.BR
// (R_t column); already at quarterly frequency, so -- unlike the baseline
// panel's coefficients -- none of these need dividing by 4.
ARRUS   =  0.3191;             // coefficient on world rate, R_t^us          (t=1.58)
ARRUS1  = -0.2344;             // coefficient on lagged world rate, R_{t-1}^us (t=-1.13)
ARY     = -0.0221;             // coefficient on output, y_t                 (t=-0.19)
ARY1    =  0.0209;             // coefficient on lagged output, y_{t-1}      (t=0.18)
ARIV    =  0.0016;             // coefficient on investment, iv_t            (t=0.03)
ARIV1   = -0.0264;             // coefficient on lagged investment, iv_{t-1} (t=-0.54)
ARTBY   = -0.2377;             // coefficient on trade balance/GDP, tby_t    (t=-1.75)
ARTBY1  = -0.0578;             // coefficient on lagged trade balance/GDP, tby_{t-1} (t=-0.41)
ARR1    =  0.8385;             // own lag coefficient                        (t=18.56)

// Estimated persistence of the world interest rate, eq. (6.22). BRAZIL
// ONLY (Table 6.1.BR, Rus_{t-1} row of the Rus_t column, N=128, t=41.27).
// Substantially more persistent than the baseline panel's 0.830391 --
// expected, given this AR(1) is estimated over 1996-2026 rather than
// 1994-2001, and the US rate's own persistence is well known to be
// regime-dependent (this sample spans the Greenspan/Bernanke/Yellen/
// Powell eras, the 2008 GFC, and the 2022-23 hiking cycle).
ARUSRUS1 = 0.9641;

// Persistence of the i.i.d. shock processes er, erus (0 = pure innovations)
RHOR    = 0;
RHORUS  = 0;

model;

// --- eq. (6.2): Cobb-Douglas technology ---------------------------------
exp(yy) = exp(k)^ALFA * exp(h)^(1-ALFA);

// --- eqs. (6.3), (6.4), (6.19): firm factor demands and shadow rate -----
#rd = exp(r)/(1-PSSI*(d-DBAR));                                   // R^d_t
#w  = (1-ALFA)*exp(k)^ALFA*exp(h)^(-ALFA) / (1+ETA*(rd-1)/rd);    // wage, eq (6.3)
#u  = ALFA*exp(k)^(ALFA-1)*exp(h)^(1-ALFA);                        // rental rate, eq (6.4)

// --- eq. (6.9): household budget constraint / debt accumulation --------
// NOTE: d (debt / net foreign liabilities) is kept as a LEVEL variable, not
// logged, because this calibration has D_ss < 0 (the economy is a net
// creditor in steady state once the working-capital wedge is netted out),
// so log(d) would be undefined. This mirrors the tb/tby treatment below.
exp(r(-1))*d(-1) + (PSSI/2)*(d-DBAR)^2 - w*exp(h) - u*exp(k)
  + exp(c) + exp(iv) = d;

// --- eq. (6.21): trade balance and TB/GDP (kept as level variables) ----
tb  = exp(yy) - exp(c) - exp(iv) - (PSSI/2)*(d-DBAR)^2;
tby = tb/exp(yy);

// --- eq. (6.7): capital accumulation with adjustment costs -------------
#ac = exp(s3)/exp(k) - (PHI/2)*(exp(s3)/exp(k)-DELTA)^2;   // Phi(s3/k)
exp(k(+1)) = (1-DELTA)*exp(k) + exp(k)*ac;

// --- eq. (6.6): investment-project pipeline -----------------------------
exp(s1) = exp(s0(-1));
exp(s2) = exp(s1(-1));
exp(s3) = exp(s2(-1));

// --- eq. (6.5): investment aggregator ------------------------------------
exp(iv) = (exp(s0)+exp(s1)+exp(s2)+exp(s3))/4;

// --- eq. (6.11): consumption Euler equation ------------------------------
exp(la(+1)) = (exp(c(+1)) - MU*exp(c) - exp(h(+1))^OMEGA/OMEGA)^(-GAMA);

// --- eq. (6.12): labor supply ---------------------------------------------
#rd_lead = exp(r(+1))/(1-PSSI*(d(+1)-DBAR));
#w_lead  = (1-ALFA)*exp(k(+1))^ALFA*exp(h(+1))^(-ALFA) / (1+ETA*(rd_lead-1)/rd_lead);
(exp(c(+1))-MU*exp(c)-exp(h(+1))^OMEGA/OMEGA)^(-GAMA)*exp(h(+1))^(OMEGA-1)
  = exp(la(+1))*w_lead;

// --- eq. (6.13): bond Euler equation ---------------------------------------
exp(la)*(1-PSSI*(d-DBAR)) = BETTA*exp(r)*exp(la(+1));

// --- eqs. (6.14)-(6.16): investment-project valuation ----------------------
exp(nu0(+1)) = 0.25;
BETTA*exp(la(+1))*exp(nu1(+1)) = (BETTA/4)*exp(la(+1)) + exp(la)*exp(nu0);
BETTA*exp(la(+1))*exp(nu2(+1)) = (BETTA/4)*exp(la(+1)) + exp(la)*exp(nu1);

// --- eq. (6.17): Tobin's Q / stage-3 project pricing ------------------------
#acp_prime_lead = 1 - PHI*(exp(s3(+1))/exp(k(+1)) - DELTA);   // Phi'(s3'/k')
BETTA*exp(la(+1))*exp(qq(+1))*acp_prime_lead
  = (BETTA/4)*exp(la(+1)) + exp(la)*exp(nu2);

// --- eq. (6.18): capital pricing equation -----------------------------------
#x3_lead    = exp(s3(+1))/exp(k(+1));
#acp_lead   = x3_lead - (PHI/2)*(x3_lead-DELTA)^2;      // Phi(s3'/k')
#acpp_lead  = 1 - PHI*(x3_lead-DELTA);                  // Phi'(s3'/k')
#u_lead     = ALFA*exp(k(+1))^(ALFA-1)*exp(h(+1))^(1-ALFA);
exp(la)*exp(qq) = BETTA*( exp(la(+1))*exp(qq(+1))*
                  (1-DELTA+acp_lead-x3_lead*acpp_lead) + exp(la(+1))*u_lead );

// --- eq. (6.22): world interest rate law of motion (calibrated, not estimated) --
rus - log(RUSSS) = ARUSRUS1*(rus(-1)-log(RUSSS)) + erus;

// --- eq. (6.20): country interest rate law of motion (calibrated, not estimated) --
r - log(RSS) = ARRUS*(rus-log(RUSSS))     + ARRUS1*(rus(-1)-log(RUSSS))
             + ARY*(yy-log(YYSS))          + ARY1*(yy(-1)-log(YYSS))
             + ARIV*(iv-log(IVSS))         + ARIV1*(iv(-1)-log(IVSS))
             + ARTBY*(tby-TBYSS)           + ARTBY1*(tby(-1)-TBYSS)
             + ARR1*(r(-1)-log(RSS))       + er;

// --- shock innovation processes --------------------------------------------
er   = RHOR*er(-1)     + eps_r;
erus = RHORUS*erus(-1) + eps_rus;

end;

// ------------------------------------------------------------------
// Steady state, computed in closed form exactly as in the model's own
// steady-state formulas (no numerical solver needed).
// ------------------------------------------------------------------
steady_state_model;
  NU0 = 1/4;
  NU1 = 1/4 + NU0/BETTA;
  NU2 = 1/4 + NU1/BETTA;
  QQ  = 1/4 + NU2/BETTA;                              // Tobin's Q, steady state

  Uss    = QQ*(RSS-1+DELTA);                          // rental rate, eq (6.18) at SS
  hoverk = (Uss/ALFA)^(1/(1-ALFA));                    // eq (6.4) at SS
  H = ( (1-ALFA)/(1+ETA*(RSS-1)/RSS)
        * (ALFA/(QQ*(RSS-1+DELTA)))^(ALFA/(1-ALFA)) )^(1/(OMEGA-1));  // (6.3)+(6.12) at SS
  K  = H/hoverk;
  IVSS = DELTA*K;                                      // eq (6.7) at SS
  YYSS = K^ALFA*H^(1-ALFA);                             // eq (6.2)
  TBss = STB*YYSS;                                      // calibration target
  Css  = YYSS - IVSS - TBss;                            // eq (6.21) at SS
  Wss  = (1-ALFA)*YYSS/H/(1+ETA*(RSS-1)/RSS);           // eq (6.3)
  Dss  = (Css+IVSS-Wss*H-Uss*K)/(1-RSS);                // eq (6.9) at SS
  TBYSS = TBss/YYSS;
  DBAR  = Dss;
  LAss  = (Css*(1-MU) - H^OMEGA/OMEGA)^(-GAMA);         // eq (6.11) at SS

  c  = log(Css);
  h  = log(H);
  k  = log(K);
  d  = Dss;
  s0 = log(IVSS); s1 = log(IVSS); s2 = log(IVSS); s3 = log(IVSS);
  iv = log(IVSS);
  yy = log(YYSS);
  tby = TBYSS;
  tb  = TBss;
  r   = log(RSS);
  rus = log(RUSSS);
  la  = log(LAss);
  qq  = log(QQ);
  nu0 = log(NU0); nu1 = log(NU1); nu2 = log(NU2);
  er  = 0;
  erus = 0;
end;

steady;
check;


// ------------------------------------------------------------------
// Shocks, for the IRF-matching step specifically. Per Uribe and Yue
// (2006, JIE), Section 3.1.4: Figure 2 ("a unit innovation in the
// country spread shock") and Figure 3 ("a one percentage point increase
// in the US interest rate shock") BOTH use a fixed shock size of 0.01 in
// raw log units -- confirmed numerically: log(1.01+0.01)-log(1.01) =
// 0.0099 =~ 0.01, exactly the "1 percentage point on a x100-scaled axis"
// convention both figures' interest-rate panels use. This is a
// deliberate, FIXED normalization for the estimation exercise itself,
// NOT the shocks' true estimated standard deviations. Both sides of the
// comparison (the model's IRFs here, and the target IRFs loaded below)
// must use this same 0.01 size for the matching to be apples-to-apples.
//
// The shocks' TRUE estimated standard deviations -- stated directly in
// the paper's text, eqs. (13) and (15) -- are 0.031 for epsilon^r and
// 0.007 for epsilon^rus. These are the right values to use for anything
// OTHER than this specific IRF-matching step (e.g. a stochastic
// simulation meant to reflect actual shock volatility, or a variance
// decomposition); they are not used here.
// ------------------------------------------------------------------
shocks;
  var eps_r;   stderr 0.01;
  var eps_rus; stderr 0.01;
end;

// Load the target IRFs produced by uribe_yue_VAR_Brazil_target_export.m.
load IRF_target_uribe_yue_Brazil;

// ------------------------------------------------------------------
// Parameters to estimate by IRF matching: the four "deep" frictions the
// paper itself calibrates against the estimated VAR's impulse responses.
// Priors are wide, uniform ranges bracketing the authors' own point
// estimates -- UNCHANGED from the baseline model.
//
// STARTING VALUES, however, are NOT the paper's own point here (unlike
// the baseline). Brazil's estimated R/Rus persistence (ARR1=0.8385,
// ARUSRUS1=0.9641) makes the model Blanchard-Kahn UNSTABLE at the paper's
// calibration. That scan found stability is overwhelmingly driven by
// MU (habit formation), with a "sweet spot" band of MU in roughly [0.47,
// 0.86]. Economically this
// is a sensible mechanism: habit formation dampens the consumption
// response to shocks, and Brazil's much more persistent estimated rate
// process apparently needs substantially more of that dampening than the
// original 7-country panel's calibration provides.
//
// The values below are the closest-to-the-paper STABLE point found by
// that scan 
// ------------------------------------------------------------------
estimated_params;
  PSSI, uniform_pdf, , , 0.000001, 0.01;    // debt-adjustment cost intensity  (paper: 0.000425)
  PHI,  uniform_pdf, , , 1,        200;     // capital-adjustment cost intensity (paper: 72.8268)
  ETA,  uniform_pdf, , , 0,        3;       // working-capital wedge share      (paper: 1.2023)
  MU,   uniform_pdf, , , 0,        0.95;    // external habit intensity         (paper: 0.2037)
end;

estimated_params_init;
  PSSI, 0.000358;
  PHI,  45.5303;
  ETA,  1.0010;
  MU,   0.4697;
end;

// ------------------------------------------------------------------
// Target IRFs: 4 variables (output, investment, trade balance/GDP, country
// rate -- NOT the world rate, which is a regressor/shock source in the
// VAR but never itself a matched target) x 2 shocks x 24 quarters = 192
// points, exactly matching the paper's own eq. (17) setup ("we are
// setting 4 parameter values to match 192 points").
//
// WEIGHTING: eq. (17) weights each of the 192 points by the inverse of
// its delta-method-estimated variance (down-weighting noisier, wide-
// error-band points). That weighting matrix is NOT implemented here --
// uribe_yue_VAR_baseline.mlx's own Section 3b notes this delta-method
// Jacobian is a substantially bigger undertaking than the point estimates
// and suggests a residual/pairs bootstrap as a pragmatic alternative;
// that bootstrap has not been built yet either. Every point below is
// therefore implicitly EQUAL-WEIGHTED, a known simplification relative to
// the paper's actual estimator. This mainly matters for how much weight
// long-horizon / noisily-estimated points get relative to sharply-
// estimated short-horizon ones -- worth revisiting if the estimated
// parameters here don't line up well with Table 3's reported values.
// ------------------------------------------------------------------
matched_irfs;
  var yy;  varexo eps_r;   periods 1:24; values(yy_eps_r);
  var iv;  varexo eps_r;   periods 1:24; values(iv_eps_r);
  var tby; varexo eps_r;   periods 1:24; values(tby_eps_r);
  var r;   varexo eps_r;   periods 1:24; values(r_eps_r);

  var yy;  varexo eps_rus; periods 1:24; values(yy_eps_rus);
  var iv;  varexo eps_rus; periods 1:24; values(iv_eps_rus);
  var tby; varexo eps_rus; periods 1:24; values(tby_eps_rus);
  var r;   varexo eps_rus; periods 1:24; values(r_eps_rus);
end;

// Method of moments estimation with IRF matching
// .
method_of_moments(
  mom_method = IRF_MATCHING,

  datafile = 'IRF_target_uribe_yue_Brazil.mat',
  order = 1,
  mode_compute = 4,
  mh_nblocks = 2,
  mh_replic = 200,
  mh_jscale = 10,
  mh_drop = 0.35,
  posterior_sampling_method = 'slice',
  plot_priors = 1
);

// ------------------------------------------------------------------
// Post-estimation: recompute IRFs at the ESTIMATED parameter values (Dynare
// writes the point estimate back into M_.params after method_of_moments).
// Still at the 0.01 shock size used for matching (shocks; block above was
// not changed), so oo_.irfs is directly comparable, point for point, to
// the target IRFs in IRF_target_uribe_yue.mat -- ready for the comparison
// plots in replicate_IRF_uribe_yue_matching.m. (irf=25 here vs. 24 in the
// matching just shows one extra quarter of decay past the matched horizon;
// it doesn't affect the estimation itself.)
// ------------------------------------------------------------------
stoch_simul(order=1, irf=25) yy tby r rus iv c h k d;
