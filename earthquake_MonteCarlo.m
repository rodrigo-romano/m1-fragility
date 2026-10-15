function dataout = earthquake_MonteCarlo(eq_stress_data,rp_bins,Nsim,weibull,baseline_stress)
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
    
    K = weibull.K;
    beta = weibull.beta;
    n =weibull.n;
    
    nYT = 50*Nsim; %total number of years.

    % log of return period at center of each bin
    for k = 1:length(rp_bins)-1
        logrp(k) = (mean(rp_bins(k:(k+1))));
    end
    rp = 10.^(logrp);
    nY = diff(10.^rp_bins);  %number of years between eq on average;
%     numE = round(nYT./nY);  % HOW many of THIS earthquake in the total period
    numE = round(nYT./rp);  % how MANY of THIS earthquake in the total period

    stress_mean = interp1(eq_stress_data.logrp,eq_stress_data.mean,logrp);
    stress_std = interp1(eq_stress_data.logrp,eq_stress_data.std,logrp);
    
    %% Setup random earthquake levels and times
    earthquake_time = zeros(sum(numE),1);
    eqrp = zeros(sum(numE),1);

    is = 1;
    for k = 1:length(nY)
        nE = numE(k); %number of EQ in total
        earthquake_time(is:(is+nE-1)) = rand(nE,1)*nYT; %WHEN it occured 
        eqrp(is:(is+nE-1)) = rp(k)*ones(nE,1); % Return period
         eqid(is:(is+nE-1)) =k;
        is = is+nE;
    end
    
    % sort by the time they occur!
    [earthquake_time, ids] = sort(earthquake_time); 
    eqrp = eqrp(ids);
    eqid = eqid(ids);
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
    if(0)
        
        F_check = 1-exp(-(K*sigmaC.^n.*teT_baseline).^beta)
    
    end
    
    teT = teT_baseline; % SHOULD START WITH 50 years at 250 psi...
    Nyears = 50;
    Te = zeros(1,Nsim);
    for k =1:length(earthquake_time)
        if(earthquake_time(k)<(kk*Nyears))
            stress_eq = (stress_mean(eqid(k))+randn(1)*stress_std(eqid(k)));
            teT =teT+ (stress_eq/sigmaC).^n.*t1_50year;  %TT/(stC)^n SUM ((st)^n)
        else
    %         Ft(kk) = 1-exp(-(K*sigmaC.^n.*teT).^beta);
            Ft(kk) = 1-exp(-(K*sigmaC.^n.*teT).^beta);
            Te(kk) = teT;
            kk =kk+1;
            teT=teT_baseline;
        end
    end
    Ft(kk) = 1-exp(-(K*sigmaC.^n.*teT).^beta);
    Te(kk) = teT;
    

    
    dataout.Te = Te;
    dataout.pf = Ft;
    dataout.pfT = mean(Ft);
    Se1 =(((-log(1-mean(Ft))).^(1/beta))/(K*50*365*24*3600)).^(1/n);
    
    dataout.Se1 = Se1; % equivalent 50 year stress from mean Ft
end