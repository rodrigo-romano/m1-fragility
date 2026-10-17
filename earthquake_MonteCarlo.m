function dataout = earthquake_MonteCarlo(eq_stress_data,rp_bins,Nsim,weibull,baseline_stress,rseed,check_mode)
% runs a brute force monte carlo using equivalent stress statistics as a function of
% return period.  These are then converted to equivlent time. Earthquakes occur at random times during many century
% span. Then for each 50 years span, equivalent times are summed up.
% eq_stress_data - struct with
%   - mean (psi)
%   - std (psi) (Guassian for now)
%   - ref_period (reference period in years) (unused)
%   - logrp (log10 of return period in years)
% rp_bins
%   - edges of bins to use in (log10 return period in years)
% Nsim
%   - number of 50 year time periods to simulate.
% weibull
%  - weibull statistics to use
%  - beta,etc
% baseline_stress (PSI) constant level over time span (includes safety
% factor)
% dataout
%  - Te (equivalent time for each 50 years span)
%  - pf ( probability of failure for each 50 years span)
%  - totalpf (total probabiltiy of failure)
%  - Se (equivalent 50 year stress)

if nargin < 7 || isempty(check_mode), check_mode = 'auto'; end
check_mode = lower(check_mode);
assert(ismember(check_mode,{'auto','always','off'}),...
    'check_mode must be auto, always or off');

use_vectorized = true;
rng(rseed);

K = weibull.K;
beta = weibull.beta;
n = weibull.n;

nYT = 50*Nsim; %total number of years.
logrp = zeros(length(rp_bins)-1,1);
% Log of return period at center of each bin
for k = 1:length(logrp)
    logrp(k) = (mean(rp_bins(k:(k+1))));
end
rp = 10.^(logrp);

if(use_vectorized == 1)
    lam  = 10.^(-rp_bins);  % Exceedance rate at bin edges
    numE = poisscount(nYT * -diff(lam));
else
    numE = round(nYT./rp);
end


%% Setup random earthquake levels and times
earthquake_time = zeros(sum(numE),1);
eqrp = zeros(size(earthquake_time));
eqid = zeros(size(earthquake_time));

is = 1;
for k = 1:length(logrp)
    nE = numE(k); %number of EQ in total
    earthquake_time(is:(is+nE-1)) = rand(nE,1)*nYT; %WHEN it occured
    eqrp(is:(is+nE-1)) = rp(k)*ones(nE,1); % Return period
    eqid(is:(is+nE-1)) = k;
    is = is+nE;
end

% Sort by the time they occur!
% Sort doesn't matter when the vectorized lifetime assessment is
% emlpoyed
[earthquake_time, ids] = sort(earthquake_time);
eqid = eqid(ids);

%% Mean and std interpolation for selected RP
stress_mean = interp1(eq_stress_data.logrp,eq_stress_data.mean,logrp);
stress_std = interp1(eq_stress_data.logrp,eq_stress_data.std,logrp);
assert(~any(isnan(stress_mean)) && ~any(isnan(stress_std)),...
    'interp1 out of range')

%% Sum up equivalent time in each period!
Ft = zeros(Nsim,1); % failure for each year?
kk = 1;

sigmaC = 250;% doesnt really matter!
t1_50year = 3600*24*365*50;

base_stress = baseline_stress; %710 gives 1% pof
teT_baseline = (base_stress/sigmaC).^n.*t1_50year;
% APPROXIMATELY 4000 PSI>.. ok..
%not 700?
% check calculation
if(false)
    F_check = 1-exp(-(K*sigmaC.^n.*teT_baseline).^beta)  %#ok<UNRCH>
end

teT = teT_baseline; % SHOULD START WITH 50 years base_stress!
Nyears = 50;
Te = zeros(1,Nsim);
cumk_per50y = zeros(size(Ft));

if(use_vectorized == 1)
    % Lifetime index of each event (1..Nsim)
    life = floor(earthquake_time(:)/Nyears) + 1;
    life = min(life, Nsim);     % guard: an event at exactly Nsim*Nyears

    % One stress draw per event, vectorized
    stress_eq = stress_mean(eqid(:)) + randn(numel(eqid),1).*stress_std(eqid(:));
    stress_eq = max(stress_eq, 0);      % A Gaussian can lead to misleading negative equiv. stress

    % Equivalent time of each event, then sum per lifetime
    dT = (stress_eq/sigmaC).^n * t1_50year;
    T  = teT_baseline + accumarray(life, dT, [Nsim 1]);   % empty lifetimes -> baseline

    Ft = 1 - exp(-(K*sigmaC.^n.*T).^beta);
    Te = T;
    k_per50y = accumarray(life, 1, [Nsim 1]);   % events per lifetime
else
    for k = 1:length(earthquake_time)
        stress_eq = (stress_mean(eqid(k))+ randn(1)*stress_std(eqid(k)));
        if(earthquake_time(k) > (kk*Nyears))
            Ft(kk) = 1-exp(-(K*sigmaC.^n.*teT).^beta);
            Te(kk) = teT;
            kk = kk+1;

            teT = teT_baseline;
            cumk_per50y(kk-1) = k;
        end
        teT = teT+ (stress_eq/sigmaC).^n.*t1_50year;
    end

    Ft(kk) = 1-exp(-(K*sigmaC.^n.*teT).^beta);
    Te(kk) = teT;
    cumk_per50y(kk) = k;
    k_per50y = diff([0;cumk_per50y]);
end

dataout.Te = Te;
dataout.pf = Ft;
dataout.pfT = mean(Ft);
dataout.k_per50y = k_per50y;
dataout.pf_se = std(Ft)/sqrt(Nsim);

Se1 =(((-log(1-mean(Ft))).^(1/beta))/(K*t1_50year)).^(1/n);
dataout.Se1 = Se1; % equivalent 50 year stress from mean Ft

dataout.pf_se = std(Ft)/sqrt(Nsim);

if ~strcmp(check_mode,'off')
    mu_exp = nYT * -diff(lam);                          % expected events per bin
    Fb     = 1-exp(-(K*sigmaC.^n.*teT_baseline).^beta); % baseline-only pf
    dataout.checks = mc_checks(numE, mu_exp, earthquake_time, ...
        dataout.k_per50y, Ft, Fb, Nsim, dataout.pf_se, dataout.pfT, check_mode);
end

end

% Poisson(mu) draws using only rand. Exact, works for any mean size.
function n = poisscount(mu)
n = zeros(size(mu));
for i = 1:numel(mu)
    m = ceil(mu(i) + 8*sqrt(mu(i)) + 20);   % oversample so the sum passes mu
    t = cumsum(-log(rand(m,1)));
    n(i) = sum(t <= mu(i));
end
end


% Sanity checks of the event generation and lifetime bookkeeping (C1, C2).
function chk = mc_checks(numE, mu_exp, t, k, Ft, Fb, Nsim, pf_se, pfT, mode)
%
% What each check tests:
%
% - Total events, events per lifetime and empty fraction compare the run to
% the expected Poisson values within 4 standard deviations. At 4σ a false
% alarm is very unlikely, so a failure points to a real problem.
% - Event conservation requires sum(k) to equal the number of generated
% events. The old loop would lose events here.
% - Empty lifetimes at baseline requires every empty lifetime to have
% exactly the baseline pf, with no leftover zeros.
% - Non-finite pf catches NaNs, for example from interp1 going out of range
% when you extend the bins in C3.
% - Relative standard error flags a run as unconverged above 20%. Change
% se_max if you want a tighter limit.
nsig   = 4;       % tolerance in standard deviations
se_max = 0.20;    % max admissible relative standard error of pfT

k   = k(:);  Ft = Ft(:);
emp = (k == 0);                       % empty lifetimes
mu_tot = sum(mu_exp(:));
lamL   = mu_tot/Nsim;                 % expected events per lifetime
p0     = exp(-lamL);                  % expected fraction of empty lifetimes
if any(emp)
    at_base = mean(abs(Ft(emp)-Fb) <= 1e-9*Fb);   % fraction exactly at baseline
else
    at_base = 1;
end

%        name                         value              expected     tolerance
rows = { 'Total events',              sum(numE),         mu_tot,      nsig*sqrt(mu_tot);
    'Events per lifetime',            mean(k),           lamL,        nsig*sqrt(lamL/Nsim);
    'Empty-lifetime fraction',        mean(emp),         p0,          nsig*sqrt(p0*(1-p0)/Nsim);
    'Event conservation (sum k)',     sum(k),            numel(t),    0;
    'Empty lifetimes at baseline pf', at_base,           1,           0;
    'Non-finite pf values',           nnz(~isfinite(Ft)),0,           0;
    'Rel. std. error of pfT',         pf_se/pfT,         0,           se_max };

ok = false(size(rows,1),1);
for i = 1:size(rows,1)
    ok(i) = abs(rows{i,2} - rows{i,3}) <= rows{i,4};
end

chk.names = rows(:,1);  chk.value = cell2mat(rows(:,2));
chk.expected = cell2mat(rows(:,3));  chk.tol = cell2mat(rows(:,4));
chk.ok = ok;  chk.pass = all(ok);

if strcmp(mode,'always') || ~chk.pass
    fprintf('--- Monte Carlo checks (Nsim = %d) ---\n', Nsim);
    for i = 1:numel(ok)
        fprintf('%-32s %12.5g   expected %10.5g +/- %-9.3g %s\n', rows{i,1}, ...
            rows{i,2}, rows{i,3}, rows{i,4}, ternary(ok(i),'PASS','FAIL'));
    end
    fprintf('pfT = %.4g +/- %.2g\n', pfT, pf_se);
end
if ~chk.pass
    warning('earthquake_MonteCarlo:checks', 'Monte Carlo checks failed: %s', ...
        strjoin(chk.names(~ok), '; '));
end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end