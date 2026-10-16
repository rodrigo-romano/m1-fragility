%
% plot_OLE_spectra
%

%% Load AECOM data
%
if(~exist('data1ole','var'))
    folder = '/Users/rromano/Workspace/gmt-data/Seismic-data/fromAECOM/Deliverables/';
    fn = 'OLE-0.1sec-Horizontal-Vertical-ScaledGMs.xlsx';
    filename = [folder,fn];
    data1ole = xlsread(filename,'A4:BL45713');
    % data1 = xlsread(filename,'A4:BL20000');
end

if(~exist('data2ole','var'))
    folder = '/Users/rromano/Workspace/gmt-data/Seismic-data/fromAECOM/Deliverables/';
    fn = '01-sec CMS_CS_UHS plots_Vs30-1500.xlsx';
    filename = [folder,fn];
    % [~,sheet]=xlsfinfo(filename);
    
    data2ole = xlsread(filename,'0.1-sec CMS &amp; CS','A4:G25');
    % 7col period (sec)	| CMS| CMS+1? | CMS+2? | CMS-1? | CMS-2? | UHS
end
Tole = data2ole(:,1);

if(~exist('data2sle','var'))
    folder = '/Users/rromano/Workspace/gmt-data/Seismic-data/fromAECOM/Deliverables/';
    fn = '01-sec CMS_CS_UHS plots_Vs30-1500.xlsx';
    filename = [folder,fn];
    % [~,sheet]=xlsfinfo(filename);
    
    data2sle = xlsread(filename,'0.1-sec CMS &amp; CS','I4:O25');
    % 7col period (sec)	| CMS| CMS+1? | CMS+2? | CMS-1? | CMS-2? | UHS
end
Tsle = data2sle(:,1);

%% Calculate OLE short-period CS
% Path to files of https://github.com/GMTO/general_IM_tools or
% https://github.com/rodrigo-romano/general_IM_tools
gen_im_tools_folder = '/Users/rromano/Workspace/fromDaveS/general_im_tools/seismic';
addpath(gen_im_tools_folder);

f = logspace(-1,2,100);
f = sort([f,10]);
damp = .05;
Sa = zeros(length(f),1,16);
g = 1;
N = 10; % Number of samples to take median over for RotD50
for k = 1:16
   ind = (k-1)*4+((1:3));
   dt = diff(data1ole(1:2,ind(1)));
   Sa(:,1,k) = SpectralPSA01(data1ole(:,ind(2:3)),f,dt,damp,g,N);
end

rmpath(gen_im_tools_folder);

%%
%
figure(40);
set(gcf,'Units','pixels','Position',[224,125,809,420])

subplot(1,2,1)
loglog(Tole,data2ole(:,7),'.-','LineWidth',1.45); hold on; 
loglog(Tole,data2ole(:,2),'.-','LineWidth',1.45);
loglog(Tsle,data2sle(:,7),'.-','LineWidth',1.45);
loglog(Tsle,data2sle(:,2),'.-','LineWidth',1.45);
xlabel('Period (s)');
ylabel('Spectral acceleration (g)');
legend('Horizontal OLE UHS','Short-period (T=0.1s) OLE CMS',...
    'Horizontal SLE UHS','Short-period (T=0.1s) SLE CMS',...
    'Location','Southwest','FontSize',10)
grid on;
set(gca,'FontSize',10);
hold off;

subplot(1,2,2)
h_uhs = loglog(Tole,data2ole(:,7),'-','LineWidth',2); hold on; 
h_cms = loglog(Tole,data2ole(:,2),'Color','#D95319','LineWidth',2);
% h_psig = loglog(T,data2ole(:,3),'r:','LineWidth',2);  
% h_msig = loglog(T,data2ole(:,5),'r:','LineWidth',2); 
h_p2sig = loglog(Tole,data2ole(:,4),'-.','Color','#D95319','LineWidth',1.45);  
h_m2sig = loglog(Tole,data2ole(:,6),'-.','Color','#D95319','LineWidth',1.45); 
h_cs = loglog(1./f, squeeze(Sa(:,1,1:16)),'Color',[.5 .5 .5],'LineWidth',0.5);

legend([h_uhs, h_cms, h_p2sig, h_cs(1)],...
    'Horizontal OLE UHS','Short-period (T=0.1s) OLE CMS',...
    'Short-period OLE CMS \pm 2\sigma','16 short-period OLE CS',...
    'Location','Southwest','FontSize',10)
xlabel('Period (s)');
ylabel('Spectral acceleration (g)');
grid on;
set(gca,'FontSize',10);
hold off;

