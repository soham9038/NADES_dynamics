%% Script to calculate Theoretical vhc %%

% Create axes
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',28,'LineWidth',3,'TickLength',[0.015 0.025]);
% set(gca,'yMinorTick','Off');  
% axes1.XScale='log';
 %axes1.YScale='log';
% xlabel('Distance','FontSize',30,'Interpreter','latex');
% y=ylabel('$P$','FontSize',42,'Interpreter','latex',...
%     'Rotation',90);
xlabel('$\Delta r (\AA)$','FontSize',36,'Interpreter','latex');
% y=ylabel('$P(W^*)$','FontSize',36,'Interpreter','latex',...
%     'Rotation',90);
y=ylabel('$4 \pi r^2 G_s(\Delta r, \tau_{NG})$','FontSize',36,'Interpreter','latex',...
    'Rotation',90);
hold on;
%xlim([0 8.0]);
%ylim([0 0.7]);
set(gca,'XMinorTick','on','YMinorTick','on');
pbaspect([1.2 1. 1.]);
format long;


inpdat=importdata('VHC_93.xyz');
pfc=inpdat;
h1=histogram(pfc(:,1),'Normalization','pdf','DisplayName','$1M-MWC$','FaceColor','m');
h1.NumBins=200; 
h1.FaceAlpha=0.3;
h1.EdgeAlpha=0.1;
pd5=fitdist(pfc(:,1),'normal');
sk1=skewness(pfc(:,1));
ku1=kurtosis(pfc(:,1));
h1.Visible='off';
h1.HandleVisibility='off';
hold on;

ydat=h1.Values;
%y1dat=ydat*4*3.14.*xdat.*xdat;
dummy=size(ydat);
len=dummy(2);
xdat=zeros(1,len);
xraw=h1.BinEdges;

for i=1:(len)
    xdat(i)=xraw(i)+((xraw(i+1)-xraw(i))/2.);
end

pl4=plot(xdat,ydat, '-sm','DisplayName','$1M-MWC$');
pl4.LineWidth=2;
pl4.MarkerSize=6;
%%%Fitting With Custom Gaussian Expression 
ft = fittype( '4*3.1415926535897932384626433832795*xdat*xdat*(a*(exp(-(3*xdat*xdat)/(2*0.0401978))))', 'independent', 'xdat', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 40.9363061498;
fp = fit(xdat(:), ydat(:), ft, opts);
pfit=plot(fp, 'c');
pfit.LineWidth=3;
pfit.LineStyle='--';
%%% Legend Updated
legend('$G_s^{Sim}(\Delta r, \tau_{NG})$', '$G_s^{Theo}(\Delta r, \tau_{NG})$');
xlabel('$\Delta r \, (\AA)$');
y=ylabel('$4 \pi r^2 G_s(\Delta r, \tau_{NG})$','FontSize',30,'Interpreter','latex',...
    'Rotation',90);
a4 = trapz(xdat,ydat);


lgd=legend('show');      
hold on;
lgd.Box='off';
lgd.FontSize=26;
lgd.Location='best';
lgd.Interpreter='latex';

fileID = fopen('r.txt','w');
fprintf(fileID,'%6.2f %12.8f\n',xdat);
fclose(fileID);

fileID = fopen('simu.txt','w');
fprintf(fileID,'%6.2f %12.8f\n',ydat);
fclose(fileID);

fileID = fopen('theo.txt','w');
fprintf(fileID,'%6.2f %12.8f\n',fp(xdat));
fclose(fileID);