% Part of worked example
%   Takes stress results from S91 and performs an example Monte Carlo
%   Study. This is meant only as an example not an exhaustive study.

% Fork from post_process_s17*.m

% Note/Issue: 16 CS earthquakes have Sa assumed to be log distirbutded. But
% this uses normal standard deviation. (doesnt look like ANY distribution
% would be a great fit)

% ISSUE:
% need better traceability of K/beta used in F computation

% NOTE:
% the "tet" (total equivalent time) variable is overwritten. In both cases
% its a 3DOF matrix. But the order of the indexes is different!

% Issue:
% Extrapolation appears to be log based when linear would be slighlty more
% accurate.

% STRESS GROUPS for reference
% 1,Faceplate ((F)Fibrefrax)
% 2,Web ends (F)
% 3,Fillet, back plate (F)
% 4,Cross-pin holes (QF150)
% 5,Backplate near venthole (F)
% 6,Backplate away form venthole (originally ALL else) ( polished except GMT1)

groupnames = {'faceplate','web ends','backplate: fillet','cross-pin holes','backplate: vent hole','GMT1 back'};
grouplimits = [714,714,714,472,714,427];
groupSF = {'Fibrefrax','Fibrefrax','Fibrefrax','QF-150','Fibrefrax','Polished, crack under loadspreader'}; 

% S91 runs for reference
% 1-16 SLE
% 17 -32 OLE
% 33-37 SLE CMS (35=1) Scale = (.54,.735,1,1.36,1.85)
% 38 -42 OLE CMS (40=1)Scale = (.54,.735,1,1.36,1.85)

% close all;

%% PART 1
% 1. Loads all 42 stress cases into "stressdata" for later use
% 2. Temporary calculation of 50 yr equivalent stress for each
% (earthquake,stress group, mirror)

Nseg =  7;
Se = 300; %normalizing stress psi for equivlanet time (ends up not having any effect)
tn = 60*60*24*365*50; % normalizing time in seconds (50 years)
S_annealing = 100; %annealing stress
S_thermal = 100; %thermal stress
n = 30.5*ones(1,Nseg);
dt = .01;
SFM = 1;
SFTA = 1;

stressdata = cell(42,1);
tet = zeros(numel(stressdata),Nseg,6);
eqs = zeros(size(tet));

for k = 1:numel(stressdata) % ALL
    % 1 Import stress data    
    % fldr = 'T:\SeismicResults\S91\dschwartz\stress';
    % fldr = '/Users/rromano/Workspace/fromDaveS/SeismicResults/S91/dschwartz/stress';
    fldr = '/home/rromano/Workspace/seismic-tool/Results/S91/stress';
    stdy = 'S91';
    
    if k<10
        strid = sprintf('00%d',k);
    else
        strid = sprintf('0%d',k);
    end
    fn = fullfile(fldr, sprintf('GMT_Earthquakes_%s_%s.mat',stdy,strid));
    load(fn);
    % Stress1 (1x6) cellarray of stress matrices. Each matrix Nx7 (index x
    % segment)
    % params_all parameters of run (not used here)
    % ids (1x7) of time ids (not used here)

    stressdata{k} = Stress1;% store each run for later!
    psi2Nm2 = 6894.744825; % not used
    
    % Calculate equivalent time and stress for each run, for given safety
    % factors and Thermal/annealing assumption
    for k2 = 1:Nseg
        for k3 =1:6 % stress groups!
            %tet = total equivalent time for a given earthquake
            %eqs = equivalent stress for that earthquake (50 yr)
           tet(k,k2,k3) = eqtime(Stress1{k3}(:,k2)*SFM+(S_annealing+S_thermal)*SFTA,n(k2),Se,dt); 
           eqs(k,k2,k3) = eqstress(tet(k,k2,k3),Se,n(k2),tn);
        end
    end
    
end
%% PLOT histograms of equivalent stress (S1)
%%
if(1 == 0) 
    % combine ALL segments
    % sids = 1; % only do 1 segment
    sids = 1:7; % all 7 segments
    for k3 = 1:6
        figure(101)
        subplot(2,3,k3);
        histogram((eqs(1:16,sids,k3)),10); hold on;
        plot(eqs(35,sids,k3),zeros(length(sids),1),'*')
        title(sprintf('Stress Group %d',k3));
        xlabel('50 yr equivalent stress (psi)')
        legend('Conditional Spectra','Mean Spectra')
        hold off;
        
        figure(102)
        subplot(2,3,k3);
        histogram((eqs(17:32,sids,k3)),10); hold on;
        plot(eqs(40,sids,k3),zeros(length(sids),1),'*')
        title(sprintf('Stress Group %d',k3));
        xlabel('50 yr equivalent stress (psi)')
        legend('Conditional Spectra','Mean Spectra')
        hold off;
    end
    figure(101);
    sgtitle('SLE runs');
    figure(102);
    sgtitle('OLE runs');
end
%% PLOT histograms of equivalent stress (S1) - max stress across segments
%%
if(1 == 0)
    % sids = 1; % only do 1 segment
    sids = 1:7; % all 7 segments
    for k3 = 1:6
        figure(201)
        set(gcf,'Units','pixels','Position',[448,259,809,420])
        subplot(2,3,k3);
        histogram(max(eqs(1:16,sids,k3),[],2),10); hold on;
        plot(max(eqs(35,sids,k3),[],2),zeros(length(sids),1),'r*')
        title(sprintf('Stress Group %d',k3));
        xlabel('50 yr equivalent stress (psi)')
        if(k3==1), legend('Conditional Spectra','Mean Spectra'); end
        hold off;
        
        figure(202)        
        subplot(2,3,k3);
        histogram(max(eqs(17:32,sids,k3),[],2),10); hold on;
        plot(max(eqs(40,sids,k3),[],2),zeros(length(sids),1),'r*')
        title(sprintf('Stress Group %d',k3));
        xlabel('50 yr equivalent stress (psi)')
        if(k3==1), legend('Conditional Spectra','Mean Spectra'); end
        hold off;
    end
    
    figure(201);
    sgtitle('SLE runs')
    set(gcf,'Units','pixels','Position',[448,259,809,420])
    figure(202);
    set(gcf,'Units','pixels','Position',[224,125,809,420])
    sgtitle('OLE runs')
end
%%
% PART II
% condition hazard data to prepare for fragilty studies
% outputs are
% eq_stress_data and eq_stress_data_fractile
% which includes return periods for each level of earthquake studied, for
% mean spectra, and 95th percentil fracitle.
%
% MAP scales from study (33-40) to return period
% pga hazard curve from: Final composite hazard curves_no sigma red_Vs30-760 & 1500mps
% col1 = g (Sa)
% col2-8 = freq of exceedance (1/return period)
% 2 = 0
% 3 = 0.1 (SHORT)
% 4 = .15
% 5 = .3
% 6 =.75
% 7 = 2 (LONG)
% 8 = 4 
%SLE = 1/2475 frequency
%OLE = 1/100 frequency
hzd=[0.014	9.03E-01	2.02E+00	2.14E+00	1.58E+00	4.81E-01	9.26E-02	2.82E-02;...
    0.021	5.22E-01	1.35E+00	1.43E+00	9.83E-01	2.74E-01	5.50E-02	1.72E-02;...
    0.035	2.43E-01	7.55E-01	8.00E-01	5.02E-01	1.33E-01	2.81E-02	8.27E-03;...
    0.0525	1.26E-01	4.48E-01	4.75E-01	2.79E-01	7.47E-02	1.59E-02	4.03E-03;...
    0.07	7.65E-02	2.99E-01	3.17E-01	1.81E-01	4.95E-02	1.02E-02	2.21E-03;...
    0.105	3.63E-02	1.60E-01	1.72E-01	9.52E-02	2.72E-02	5.01E-03	8.13E-04;...
    0.14	2.06E-02	9.94E-02	1.07E-01	5.94E-02	1.73E-02	2.83E-03	3.55E-04;...
    0.175	1.30E-02	6.72E-02	7.32E-02	4.07E-02	1.18E-02	1.73E-03	1.74E-04;...
    0.21	8.80E-03	4.81E-02	5.28E-02	2.95E-02	8.52E-03	1.12E-03	9.23E-05;...
    0.28	4.60E-03	2.76E-02	3.09E-02	1.73E-02	4.82E-03	5.25E-04	3.10E-05;...
    0.35	2.69E-03	1.76E-02	2.00E-02	1.11E-02	2.96E-03	2.74E-04	1.23E-05;...
    0.42	1.68E-03	1.19E-02	1.37E-02	7.56E-03	1.92E-03	1.55E-04	5.48E-06;...
    0.49	1.11E-03	8.49E-03	9.91E-03	5.37E-03	1.30E-03	9.23E-05	2.67E-06;...
    0.56	7.56E-04	6.26E-03	7.40E-03	3.94E-03	9.09E-04	5.76E-05	1.39E-06;...
    0.7	3.84E-04	3.67E-03	4.44E-03	2.27E-03	4.78E-04	2.49E-05	4.39E-07;...
    1.05	9.75E-05	1.27E-03	1.60E-03	7.49E-04	1.29E-04	4.55E-06	4.29E-08;...
    1.4	3.28E-05	5.44E-04	7.13E-04	3.09E-04	4.49E-05	1.18E-06	6.10E-09;...
    1.75	1.31E-05	2.67E-04	3.59E-04	1.46E-04	1.84E-05	3.79E-07	9.54E-10;...
    2.1	5.92E-06	1.43E-04	1.97E-04	7.58E-05	8.47E-06	1.41E-07	9.17E-11;...
    2.8	1.55E-06	4.91E-05	7.02E-05	2.48E-05	2.27E-06	2.56E-08	0.00E+00];

oleSa = 10^interp1(log10(hzd(:,3)),log10(hzd(:,1)),log10(1/100));
sleSa = 10^interp1(log10(hzd(:,3)),log10(hzd(:,1)),log10(1/2475));

% use OLE for first 4 return periods
rp_study(1:4) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(oleSa*[.54,.735,1,1.36])));
% Use SLE for last 5
rp_study(5:9) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(sleSa*[.54,.735,1,1.36,1.85]),'linear','extrap'));

eq_stress_data.ref_period = 50;% years...unused
eq_stress_data.logrp = log10(rp_study);


%%
% calculate return periods. Worst case fractile ~ 30% higher than mean.
% 
eq_stress_data_fractile = eq_stress_data;
oleSa_ = oleSa/1.3;
sleSa_ = sleSa/1.3;% look on hazard curve for return period associated with 30% lower magnitude earthquake.
% use OLE for first 4 return periods
rp_study(1:4) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(oleSa_*[.54,.735,1,1.36])));
% Use SLE for last 5
rp_study(5:9) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(sleSa_*[.54,.735,1,1.36,1.85]),'linear','extrap'));
eq_stress_data_fractile.logrp = log10(rp_study);

%%
% check max over 7 mirrors vs. size of eq
figure(301);
if(~exist('sids','var')), sids = 1:7; end
loglog([.54,.735,1,1.36,1.85]*sleSa, (max(eqs(33:37,sids,k3),[],2)),'^--',...
    [.54,.735,1,1.36,1.85]*oleSa, (max(eqs(38:42,sids,k3),[],2)),'v--')
xlabel('T=0.1s Ground SA (g)')
ylabel('Max 50 yr eq. stress over 7 segments (psi)')
legend('Scaled SLE mean spectra runs','Scaled OLE mean spectra runs')
grid on;

%% plot some hazard stuff
if(0)
    figure(103); %#ok<UNRCH>
    N1 = size(hzd,1);
    loglog(hzd(:,1),hzd(:,3),'k',...
        hzd(:,1)*1.3,hzd(:,3),'r'...
        ,hzd(:,1),ones(N1,1)*1/100,'k:',hzd(:,1),ones(N1,1)*1/2475,'k:','LineWidth',2);
    grid on;
    xlabel('Spectral Accel (g) T = 0.1s');
    ylabel('Frequency of Exceedance (times/year)')
    legend('mean hazard','95 percentile hazard','OLE/SLE frequencies')
end

%% PART III
% 
% now take Stressdata{} and convert to eq time for given set of safety
% factors
clear Se1 Se2
base_mech = 50;% Assumed long term mechanical stress (applied for 50 yr)
base_thermal_annel = 200;   %psi
sf_thermal_anneal = 1;
sf_mechanical = 1:.2:5;
% sf_mechanical = 3.5;%2; % * <---
% weights = [7,3,2]/12; % used if ZA0,ZA30,ZA60

Nsim = 5000; %1e5;%
% rp_bins = 1:.1:3.4;
rp_bins = 1.7:.1:3.4; % start higher with no RLE % * <---

weibull.n = 30.5;
weibull.beta = .4522;
weibull.K = 3e-101;
clear pof1 pof2 req_lim1 req_lim2
Se1 = zeros(length(sf_thermal_anneal), length(sf_mechanical), 6);
Se1_fract = zeros(size(Se1));
for k1 = 1:length(sf_thermal_anneal)
    for k2=1:length(sf_mechanical)
        clear tet
        % collect eqtime for given safety factors
        for j = 1:42
            for k = 1:Nseg
                for k3=1:6
                    s = stressdata{j}{k3}(:,k)*sf_mechanical(k2)+...
                        sf_thermal_anneal(k1)*base_thermal_annel;
                    tet(j,k3,k) = eqtime(s,n(k),Se,dt);
                end
            end
        end
        
        % Convert to eq stress
        stress_eq = eqstress(tet,Se,weibull.n,tn);
        % (tet/(3600*24*365*50)).^(1/n(1))*Se; % use 50 year stress!        
        stress_eq = max(stress_eq,[],3); % max across segments
        
        S1 = stress_eq;
        for k3 =1:6
            fprintf('Monte Carlo run of group %s\n',groupnames{k3});
            % now build mean/sd for various rp
            % -- OLE ---
            sdratio = std(S1(16+(1:16),k3))/mean(S1((1:16)+16,k3));
            mean1 = S1(32+5+(1:4),k3);
            eq_stress_data.mean(1:4)=mean1;
            eq_stress_data.std(1:4) = sdratio*mean1;
            % -- SLE --- 
            sdratio = std(S1(1:16,k3))/mean(S1(1:16,k3));
            mean1 = S1(32+(1:5),k3);
            eq_stress_data.mean(5:9)=mean1;
            eq_stress_data.std(5:9) = sdratio*mean1;
            % MC
            dataout = earthquake_MonteCarlo(eq_stress_data,rp_bins,Nsim,weibull,...
                (sf_thermal_anneal(k1)*base_thermal_annel+sf_mechanical(k2)*base_mech));
%             pof1(k1,k2,k3)=dataout.pfT;
            Se1(k1,k2,k3) = dataout.Se1;
%             total_sf1(k1,k2,k3) = (410)/dataout.Se1;  
            
            % -- Fractile --
            eq_stress_data_fractile.mean = eq_stress_data.mean;
            eq_stress_data_fractile.std = eq_stress_data.std;
            dataout = earthquake_MonteCarlo(eq_stress_data_fractile,rp_bins,Nsim,weibull,...
                (sf_thermal_anneal(k1)*base_thermal_annel+sf_mechanical(k2)*base_mech));
%             pof1(k1,k2,k3)=dataout.pfT;
            Se1_fract(k1,k2,k3) = dataout.Se1; % fractile result
%             total_sf1(k1,k2,k3) = (410)/dataout.Se1;  

            if(abs(sf_mechanical(k2)-3.5)<.1 && k3 == 3)
               figure(108);
               set(gcf,'Units','pixels','Position',[224,125,809,360])
               subplot(121);
               histogram(log10(dataout.Te/3600/24/365))
               xlabel('Log_{10} of equivelant time (years)')
               subplot(122);
               histogram(log10(dataout.pf))
               xlabel('Log_{10} of probability of failure')               
               fprintf('Done!\n');
            end
        end
    end
end

%%
temp(:,:)=Se1(1,:,:);
temp2(:,:)=Se1_fract(1,:,:);
figure(1002); set(gcf,'Units','pixels','Position',[448,259,809,420])
for k = 1:6
    subplot(2,3,k)
    plot(sf_mechanical,temp(:,k),sf_mechanical,temp2(:,k),[min(sf_mechanical),max(sf_mechanical)],ones(1,2)*grouplimits(k));%,sf_mechanical,Se2);
    legend('mean hazard','95% fractile hazard',sprintf('limit:%s',groupSF{k}));
    grid on;
    % legend('with actuators','without actuators')
    ylabel('equivalent 50 year stress (all stress) psi')
    xlabel('mechanical safety factor')
    % title('Total 50 year equivalent stress (200psi+mech sf* (gravity+seismic))')
    title(groupnames{k});
    % ylim([200 800]);
    axis tight
end


%%
function t2 = eqtime(S2,n,S1,t1)
% equivalent time calculation
%The equivalent time, t2, will be the time by which the probability of failure at S2 
%  is equal to the probability of failure at S1 after an exposure of t1
% t1 is effectively the the sample time.  
% 
%
% Stress is stress (psi or N/m^2)
% n is Weibull n
% Se is normalizing stress (psi or N/m^2)
% t1 is normalizing time (s)
% S1 is sterss vector
    t2 = sum((S2/S1).^n*t1);
end
%%
function S1 = eqstress(tet,Se,n,tn)
% equivalent stress calculation
% (use a normalizing TIME)
% 5 minute?
%
% Stress is stress (psi or N/m^2)
% n is Weibull n
% Se is normalizing stress (psi or N/m^2)
% tn is normalizing time (s)
% tet is equivalent time
%     t2 = sum((S2/S1).^n*t1);
S1 = (tet/tn).^(1/n)*Se;
end