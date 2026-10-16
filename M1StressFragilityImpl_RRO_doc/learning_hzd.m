% PART II
% condition hazard data to prepare for fragilty studies
% outputs are
% eq_stress_data and eq_stress_data_fractile
% which includes return periods for each level of earthquake studied, for
% mean spectra, and 95th percentil fractile.
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
%RLE = 1/2 frequency
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

% use OLE for first 5 return periods
rp_study(1:5) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(oleSa*[.54,.735,1,1.36,1.85])));
% Use SLE for last 5
rp_study(6:10) = 1./(10.^interp1(log10(hzd(:,1)),log10(hzd(:,3)),log10(sleSa*[.54,.735,1,1.36,1.85]),'linear','extrap'));

eq_stress_data.ref_period = 50; % Years
eq_stress_data.logrp = log10(rp_study);

%% plot some hazard stuff
figure(102);
N1 = size(hzd,1);
h1 = loglog(hzd(:,1),hzd(:,2),'.-',...
    hzd(:,1),hzd(:,3),'.-',...
    hzd(:,1),hzd(:,4),'.-',...
    hzd(:,1),hzd(:,5),'.-',...
    hzd(:,1),hzd(:,6),'.-',...
    hzd(:,1),hzd(:,7),'.-',...
    hzd(:,1),hzd(:,8),'.-',...
    'LineWidth',1.45);
grid on;
xlabel('Spectral Acceleration (g)');
ylabel('Frequency of Exceedance (times/year)')
hold on

loglog(hzd(:,1),ones(N1,1)*1/2,'k--',...
    hzd(:,1),ones(N1,1)*1/100,'k--',...
    hzd(:,1),ones(N1,1)*1/2475,'k--');
text(mean(hzd(1:2,1)),1/2,'RLE','VerticalAlignment','Bottom')
text(mean(hzd(1:2,1)),1/100,'OLE','VerticalAlignment','Bottom')
text(mean(hzd(1:2,1)),1/2475,'SLE','VerticalAlignment','Bottom')
hold off;
legend(h1,'T \rightarrow \infty','T=0.1s','T=0.15s',...
    'T=0.3s','T=0.75s','T=2s','T=4s');
xlim([min(hzd(:,1)), 1])
ylim([1e-5, 2])


%% plot some hazard stuff
figure(103);
N1 = size(hzd,1);
h1 = loglog(hzd(:,1),hzd(:,3),'.-','Color','#D95319','LineWidth',1.5);
hold on
h1_ = loglog(hzd(:,1)*1.3,hzd(:,3),'.--','Color','#FF0000','LineWidth',1.5);
grid on;
xlabel('T=0.1s Spectral Acceleration (g)');
ylabel('Frequency of Exceedance (times/year)')
h1b = loglog(oleSa*[.54,.735,1,1.36,1.85], 1./rp_study(1:5),'o','Color',[.4,.4,.4],'Markersize',9);
h1c = loglog(sleSa*[.54,.735,1,1.36,1.85], 1./rp_study(6:10),'s','Markersize',9);
loglog([hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/2,'k--',...
    [hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/100,'k--',...
    [hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/2475,'k--');
text(mean(hzd(1:2,1)),1/2,'RLE','VerticalAlignment','Bottom','Fontsize',10)
text(mean(hzd(1:2,1)),1/100,'OLE','VerticalAlignment','Bottom','Fontsize',10)
text(mean(hzd(1:2,1)),1/2475,'SLE','VerticalAlignment','Bottom','Fontsize',10)
hold off;
legend([h1; h1_; h1b; h1c],'Mean hazard curve','1.3 Mean hazard',...
    'Scaled OLE CMS runs','Scaled SLE CMS runs','Fontsize',10);
axis tight

%% ZOOM plot
figure(104);
N1 = size(hzd,1);
h2 = semilogx(hzd(:,1),hzd(:,3),'.-',...
    hzd(:,1)*1.3,hzd(:,3),'.-',...
    [hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/2,'k--',...
    [hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/100,'k--',...
    [hzd(:,1);hzd(end,1)*1.3],ones(N1+1,1)*1/2475,'k--','LineWidth',1.45);
grid on;
xlabel('T=0.1s Spectral Acceleration (g)');
ylabel('Frequency of Exceedance (times/year)')

hold on
h2b = semilogx(oleSa*[.54,.735,1,1.36], 1./rp_study(1:4),'o');
h2c = semilogx(sleSa*[.54,.735,1,1.36,1.85], 1./rp_study(5:9),'s');
hold off;
ylim([-0.001, 0.012])
xlim([0.1, 2*sleSa])
legend([h2(1:3); h2b; h2c],'Mean hazard','95 percentile hazard','RLE/OLE/SLE return freq.',...
    'OLE study points', 'SLE study points')