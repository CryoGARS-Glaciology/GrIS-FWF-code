clearvars; close all;
iceberg_path = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-meltrates/';
figure_path = ['/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/figures/iceberg-melt/'];

%% Section 1: Load csv files in table format

cd(iceberg_path);
% DJG Files
Datatabledate1 = readtable('DJG_20190418-20190423_iceberg_meltinfo.csv');
Datatabledate2 = readtable('DJG_20200505-20200606_iceberg_meltinfo.csv');
% ISS Files
Datatabledate3 = readtable('ISS_20180418-20180508_iceberg_meltinfo.csv'); 
Datatabledate4 = readtable('ISS_20200318-20200404_iceberg_meltinfo.csv'); 
% HLG Files
Datatabledate5 = readtable('HLG_20160706-20160729_iceberg_meltinfo.csv'); 
Datatabledate6 = readtable('HLG_20210715-20210727_iceberg_meltinfo.csv');
% ASS Files
Datatabledate7 = readtable('ASS_20160427-20160603_iceberg_meltinfo.csv'); 
Datatabledate8 = readtable('ASS_20190327-20190417_iceberg_meltinfo.csv'); 
% NOG Files
Datatabledate9 = readtable('NOG_20210525-20210605_iceberg_meltinfo.csv');
% KNS Files
Datatabledate10 = readtable('KNS_20110320-20110402_iceberg_meltinfo.csv');
Datatabledate11 = readtable('KNS_20140525-20140614_iceberg_meltinfo.csv');
% SEK Files
Datatabledate12 = readtable('SEK_20110319-20110405_iceberg_meltinfo.csv');
Datatabledate13 = readtable('SEK_20130212-20130317_iceberg_meltinfo.csv');
Datatabledate14 = readtable('SEK_20130317-20130403_iceberg_meltinfo.csv');
Datatabledate15 = readtable('SEK_20140330-20140419_iceberg_meltinfo.csv');

%% Section 2:
% Grabbing Variable from CSV Files
% Variable below Grab DJG Data
VolumeChange1 = Datatabledate1{:,16};
VolumeChangeUncer1 = Datatabledate1{:,17};
Draft1 = Datatabledate1{:,18};
Draftrange1 = Datatabledate1{:,19};
SubmergedArea1 = Datatabledate1{:,22};
SubmergedAreaUncer1 = Datatabledate1{:,23};

VolumeChange2 = Datatabledate2{:,16};
VolumeChangeUncer2 = Datatabledate2{:,17};
Draft2 = Datatabledate2{:,18};
Draftrange2 = Datatabledate2{:,19};
SubmergedArea2 = Datatabledate2{:,22};
SubmergedAreaUncer2 = Datatabledate2{:,23};

% Variable below Grab ISS Data
VolumeChange3 = Datatabledate3{:,16};
VolumeChangeUncer3 = Datatabledate3{:,17};
Draft3 = Datatabledate3{:,18};
Draftrange3 = Datatabledate3{:,19};
SubmergedArea3 = Datatabledate3{:,22};
SubmergedAreaUncer3 = Datatabledate3{:,23};

VolumeChange4 = Datatabledate4{:,16};
VolumeChangeUncer4 = Datatabledate4{:,17};
Draft4 = Datatabledate4{:,18};
Draftrange4 = Datatabledate4{:,19};
SubmergedArea4 = Datatabledate4{:,22};
SubmergedAreaUncer4 = Datatabledate4{:,23};

% Variable below Grab HLG Data
VolumeChange5 = Datatabledate5{:,16};
VolumeChangeUncer5 = Datatabledate5{:,17};
Draft5 = Datatabledate5{:,18};
Draftrange5 = Datatabledate5{:,19};
SubmergedArea5 = Datatabledate5{:,22};
SubmergedAreaUncer5 = Datatabledate5{:,23};

VolumeChange6 = Datatabledate6{:,16};
VolumeChangeUncer6 = Datatabledate6{:,17};
Draft6 = Datatabledate6{:,18};
Draftrange6 = Datatabledate6{:,19};
SubmergedArea6 = Datatabledate6{:,22};
SubmergedAreaUncer6 = Datatabledate6{:,23};

% Variable below Grab ASS Data
VolumeChange7 = Datatabledate7{:,16};
VolumeChangeUncer7 = Datatabledate7{:,17};
Draft7 = Datatabledate7{:,18};
Draftrange7 = Datatabledate7{:,19};
SubmergedArea7 = Datatabledate7{:,22};
SubmergedAreaUncer7 = Datatabledate7{:,23};

VolumeChange8 = Datatabledate8{:,16};
VolumeChangeUncer8 = Datatabledate8{:,17};
Draft8 = Datatabledate8{:,18};
Draftrange8 = Datatabledate8{:,19};
SubmergedArea8 = Datatabledate8{:,22};
SubmergedAreaUncer8 = Datatabledate8{:,23};

% Variable below Grab NOG Data
VolumeChange9 = Datatabledate9{:,16};
VolumeChangeUncer9 = Datatabledate9{:,17};
Draft9 = Datatabledate9{:,18};
Draftrange9 = Datatabledate9{:,19};
SubmergedArea9 = Datatabledate9{:,22};
SubmergedAreaUncer9 = Datatabledate9{:,23};

% IN PROGRESS 
% % Variable below Grab KNS Data
% VolumeChange10 = Datatabledate10{:,16};
% VolumeChangeUncer10 = Datatabledate10{:,17};
% Draft10 = Datatabledate10{:,18};
% Draftrange10 = Datatabledate10{:,19};
% SubmergedArea10 = Datatabledate10{:,22};
% SubmergedAreaUncer10 = Datatabledate10{:,23};
% 
% VolumeChange11 = Datatable
% date11{:,16};
% VolumeChangeUncer11 = Datatabledate11{:,17};
% Draft11 = Datatabledate11{:,18};
% Draftrange11 = Datatabledate11{:,19};
% SubmergedArea11 = Datatabledate11{:,22};
% SubmergedAreaUncer11 = Datatabledate11{:,23};
% 
% % Variable below Grab SEK Data
% VolumeChange12 = Datatabledate12{:,16};
% VolumeChangeUncer12 = Datatabledate12{:,17};
% Draft12 = Datatabledate12{:,18};
% Draftrange12 = Datatabledate12{:,19};
% SubmergedArea12 = Datatabledate12{:,22};
% SubmergedAreaUncer12 = Datatabledate12{:,23};
% 
% VolumeChange13 = Datatabledate11{:,16};
% VolumeChangeUncer11 = Datatabledate11{:,17};
% Draft11 = Datatabledate11{:,18};
% Draftrange11 = Datatabledate11{:,19};
% SubmergedArea11 = Datatabledate11{:,22};
% SubmergedAreaUncer11 = Datatabledate11{:,23};
% 
% VolumeChange11 = Datatabledate11{:,16};
% VolumeChangeUncer11 = Datatabledate11{:,17};
% Draft11 = Datatabledate11{:,18};
% Draftrange11 = Datatabledate11{:,19};
% SubmergedArea11 = Datatabledate11{:,22};
% SubmergedAreaUncer11 = Datatabledate11{:,23};
% 
% VolumeChange15 = Datatabledate15{:,16};
% VolumeChangeUncer15 = Datatabledate15{:,17};
% Draft15 = Datatabledate15{:,18};
% Draftrange15 = Datatabledate15{:,19};
% SubmergedArea15 = Datatabledate15{:,22};
% SubmergedAreaUncer15 = Datatabledate15{:,23};

% Calculating curtain variables
% Units/Day m^3/day 
MeltFlux1 = (VolumeChange1./SubmergedArea1); 
MeltFlux2 = (VolumeChange2./SubmergedArea2); 
MeltFlux3 = (VolumeChange3./SubmergedArea3);
MeltFlux4 = (VolumeChange4./SubmergedArea4);
MeltFlux5 = (VolumeChange5./SubmergedArea5);
MeltFlux6 = (VolumeChange6./SubmergedArea6);
MeltFlux7 = (VolumeChange7./SubmergedArea7);
MeltFlux8 = (VolumeChange8./SubmergedArea8);
MeltFlux9 = (VolumeChange9./SubmergedArea9);
%% Section 3:Calculating Melt flux Error
% Calculation
m1 = VolumeChange1./SubmergedArea1;
SigmaA1 = (SubmergedAreaUncer1)./2;
SigmaM1 = abs(m1).* sqrt((VolumeChangeUncer1./VolumeChange1).^2+(SigmaA1./SubmergedArea1).^2);

m2 = VolumeChange2./SubmergedArea2;
SigmaA2 = (SubmergedAreaUncer2)./2;
SigmaM2 = abs(m2).* sqrt((VolumeChangeUncer2./VolumeChange2).^2+(SigmaA2./SubmergedArea2).^2);

m3 = VolumeChange3./SubmergedArea3;
SigmaA3 = (SubmergedAreaUncer3)./2;
SigmaM3 = abs(m3).* sqrt((VolumeChangeUncer3./VolumeChange3).^2+(SigmaA3./SubmergedArea3).^2);

m4 = VolumeChange4./SubmergedArea4;
SigmaA4 = (SubmergedAreaUncer4)./2;
SigmaM4 = abs(m4).* sqrt((VolumeChangeUncer4./VolumeChange4).^2+(SigmaA4./SubmergedArea4).^2);

m5 = VolumeChange5./SubmergedArea5;
SigmaA5 = (SubmergedAreaUncer5)./2;
SigmaM5 = abs(m5).* sqrt((VolumeChangeUncer5./VolumeChange5).^2+(SigmaA5./SubmergedArea5).^2);

m6 = VolumeChange6./SubmergedArea6;
SigmaA6 = (SubmergedAreaUncer6)./2;
SigmaM6 = abs(m6).* sqrt((VolumeChangeUncer6./VolumeChange6).^2+(SigmaA6./SubmergedArea6).^2);

m7 = VolumeChange7./SubmergedArea7;
SigmaA7 = (SubmergedAreaUncer7)./2;
SigmaM7 = abs(m7).* sqrt((VolumeChangeUncer7./VolumeChange7).^2+(SigmaA7./SubmergedArea7).^2);

m8 = VolumeChange8./SubmergedArea8;
SigmaA8 = (SubmergedAreaUncer8)./2;
SigmaM8 = abs(m8).* sqrt((VolumeChangeUncer8./VolumeChange8).^2+(SigmaA8./SubmergedArea8).^2);

m9 = VolumeChange9./SubmergedArea9;
SigmaA9 = (SubmergedAreaUncer9)./2;
SigmaM9 = abs(m9).* sqrt((VolumeChangeUncer9./VolumeChange9).^2+(SigmaA9./SubmergedArea9).^2);

%% Section 4: Plotting Melt rate vs. Draft graphs with error bars 
% Color: DJG black (k), ISS red (r), HLG green (g), ASS blue (b), NOG
% magenta (m),
% Grouping Notes; Number Following abrev correlate with order {ISS4 & ASS8}
% {ISS3 & DJG1} {DJG2 & ASS7 & NOG9} {HLG5 & HLG6}
% errorbar(x,y,yneg,ypos,xneg,xpos)
% errorbar(x,y,neg,pos)
% errorbar(Draft2,Delt_V2,'sr','markerfacecolor','r');
% plot(Draft2,Delt_V2,'sr','markerfacecolor','r');

figure(1), clf
set(gcf,'position',[50 50 1200 600]); 
subplot(1,4,1); set(gca,'fontsize',16);
errorbar(Draft4,m4,SigmaM4,SigmaM4,Draftrange4./2,Draftrange4./2,'sr','MarkerSize',5,'MarkerFaceColor','red'); hold on;

subplot(1,4,1); set(gca,'fontsize',16);
errorbar(Draft8,m8,SigmaM8,SigmaM8,Draftrange8./2,Draftrange8./2,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

xlabel('Draft (m)','fontsize',16)
ylabel('Melt rate (m/d)','fontsize',16)
axis tight
set(gca,'xlim',[0 400],'ylim',[0 1.25]);
legend({'20200318-20200404','20190327-20190417'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(Draft1,m1,SigmaM1,SigmaM1,Draftrange1./2,Draftrange1./2,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(Draft3,m3,SigmaM3,SigmaM3,Draftrange3./2,Draftrange3./2,'sr','MarkerSize',5,'MarkerFaceColor','red'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 1.25]);
legend({'20190418-20190423','20180418-20180508'},'Location','northwest','fontsize',14)
grid on
hold off


subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft2,m2,SigmaM2,SigmaM2,Draftrange2./2,Draftrange2./2,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft7,m7,SigmaM7,SigmaM7,Draftrange7./2,Draftrange7./2,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft9,m9,SigmaM9,SigmaM9,Draftrange9./2,Draftrange9./2,'sm','MarkerSize',5,'MarkerFaceColor','magenta'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 1.25]);
legend({'20200505-20200606','20160427-20160603', '20210525-20210605'},'Location','northwest','fontsize',14)
grid on
hold off


subplot(1,4,4); set(gca,'fontsize',16);
errorbar(Draft5,m5,SigmaM5,SigmaM5,Draftrange5./2,Draftrange5./2,'sg','MarkerSize',5,'MarkerFaceColor','green'); hold on;

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(Draft6,m6,SigmaM6,SigmaM6,Draftrange6./2,Draftrange6./2,'dg','MarkerSize',5,'MarkerFaceColor','green'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 1.25]);
legend({'20160706-20160729','20210715-20210727'},'Location','northwest','fontsize',14)
grid on
hold off

saveas(gcf,[figure_path,'MeltRatevsDraft.eps'],'epsc');
% sgtitle('Meltwater Rate')


%% Section 5: Plotting Melt flux Graphs vs. Draft with error bars
%OMITS LARGEST FLUX IN LAST PLOT

figure(2), clf
set(gcf,'position',[100 50 1200 600]); 
subplot(1,4,1); set(gca,'fontsize',16);
errorbar(Draft4,VolumeChange4./86400,VolumeChangeUncer4./86400,VolumeChangeUncer4./86400,Draftrange4./2,Draftrange4./2,'sr','MarkerSize',5,'MarkerFaceColor','red'); hold on;

subplot(1,4,1); set(gca,'fontsize',16);
errorbar(Draft8,VolumeChange8./86400,VolumeChangeUncer8./86400,VolumeChangeUncer8./86400,Draftrange8./2,Draftrange8./2,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

xlabel('Draft (m)','fontsize',16)
ylabel('Melt flux (m^3/s)','fontsize',16)
axis tight
set(gca,'xlim',[0 400],'ylim',[0 7]);
legend({'20200318-20200404','20190327-20190417'},'Location','north','fontsize',14)
grid on
hold off

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(Draft1,VolumeChange1./86400,VolumeChangeUncer1./86400,VolumeChangeUncer1./86400,Draftrange1./2,Draftrange1./2,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(Draft3,VolumeChange3./86400,VolumeChangeUncer3./86400,VolumeChangeUncer3./86400,Draftrange3./2,Draftrange3./2,'sr','MarkerSize',5,'MarkerFaceColor','red'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 7]);
legend({'20190418-20190423','20180418-20180508'},'Location','north','fontsize',14)
grid on
hold off


subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft2,VolumeChange2./86400,VolumeChangeUncer2./86400,VolumeChangeUncer2./86400,Draftrange2./2,Draftrange2./2,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft7,VolumeChange7./86400,VolumeChangeUncer7./86400,VolumeChangeUncer7./86400,Draftrange7./2,Draftrange7./2,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(Draft9,VolumeChange9./86400,VolumeChangeUncer9./86400,VolumeChangeUncer9./86400,Draftrange9./2,Draftrange9./2,'sm','MarkerSize',5,'MarkerFaceColor','magenta'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 7]);
legend({'20200505-20200606','20160427-20160603', '20210525-20210605'},'Location','north','fontsize',14)
grid on
hold off


subplot(1,4,4); set(gca,'fontsize',16);
errorbar(Draft5,VolumeChange5./86400,VolumeChangeUncer5./86400,VolumeChangeUncer5./86400,Draftrange5./2,Draftrange5./2,'sg','MarkerSize',5,'MarkerFaceColor','green'); hold on;

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(Draft6,VolumeChange6./86400,VolumeChangeUncer6./86400,VolumeChangeUncer6./86400,Draftrange6./2,Draftrange6./2,'dg','MarkerSize',5,'MarkerFaceColor','green'); 

xlabel('Draft (m)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 400],'ylim',[0 7]);
legend({'20190418-20190423','20180418-20180508'},'Location','north','fontsize',14)
grid on
hold off

saveas(gcf,[figure_path,'MeltFluxvsDraft.eps'],'epsc');
% sgtitle('Meltwater Flux')


%% Section 6: Graphs of Melt rate vs. Suberged area 
figure(3), clf
set(gcf,'position',[150 50 1200 600]); 
subplot(1,4,1); set(gca,'fontsize',16);
errorbar(SubmergedArea4./10^6,m4,SigmaM4,SigmaM4,SubmergedAreaUncer4./10^6,SubmergedAreaUncer4./10^6,'sr','MarkerSize',5,'MarkerFaceColor','red'); hold on;

subplot(1,4,1); set(gca,'fontsize',16);
errorbar(SubmergedArea8./10^6,m8,SigmaM8,SigmaM8,SubmergedAreaUncer8./10^6,SubmergedAreaUncer8./10^6,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

xlabel('Submerged area (km^2)','fontsize',16)
ylabel('Melt rate (m/d)','fontsize',16)
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 1.25],'xtick',[0:0.5:2.5]);
legend({'20200318-20200404','20190327-20190417'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(SubmergedArea1./10^6,m1,SigmaM1,SigmaM1,SubmergedAreaUncer1./10^6,SubmergedAreaUncer1./10^6,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(SubmergedArea3./10^6,m3,SigmaM3,SigmaM3,SubmergedAreaUncer3./10^6,SubmergedAreaUncer3./10^6,'sr','MarkerSize',5,'MarkerFaceColor','red'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 1.25],'xtick',[0:0.5:2.5]);
legend({'20190418-20190423','20180418-20180508'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea2./10^6,m2,SigmaM2,SigmaM2,SubmergedAreaUncer2./10^6,SubmergedAreaUncer2./10^6,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea7./10^6,m7,SigmaM7,SigmaM7,SubmergedAreaUncer7./10^6,SubmergedAreaUncer7./10^6,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea9./10^6,m9,SigmaM9,SigmaM9,SubmergedAreaUncer9./10^6,SubmergedAreaUncer9./10^6,'sm','MarkerSize',5,'MarkerFaceColor','magenta'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 1.25],'xtick',[0:0.5:2.5]);
legend({'20200505-20200606','20160427-20160603', '20210525-20210605'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(SubmergedArea5./10^6,m5,SigmaM5,SigmaM5,SubmergedAreaUncer5./10^6,SubmergedAreaUncer5./10^6,'sg','MarkerSize',5,'MarkerFaceColor','green'); hold on;

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(SubmergedArea6./10^6,m6,SigmaM6,SigmaM6,SubmergedAreaUncer6./10^6,SubmergedAreaUncer6./10^6,'dg','MarkerSize',5,'MarkerFaceColor','green'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt rate (m/d)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 1.25],'xtick',[0:0.5:2.5]);
legend({'20160706-20160729','20210715-20210727'},'Location','northwest','fontsize',14)
grid on
hold off

saveas(gcf,[figure_path,'MeltRatevsSubmergedArea.eps'],'epsc');
% sgtitle('Meltwater Rate')

%% Section 6: Melt flux vs. Submerged area
figure(4), clf
set(gcf,'position',[200 50 1200 600]); 
subplot(1,4,1); set(gca,'fontsize',16);
errorbar(SubmergedArea4./10^6,VolumeChange4./86400,VolumeChangeUncer4./86400,VolumeChangeUncer4./86400,SubmergedAreaUncer4./10^6,SubmergedAreaUncer4./10^6,'sr','MarkerSize',5,'MarkerFaceColor','red'); hold on;

subplot(1,4,1); set(gca,'fontsize',16);
errorbar(SubmergedArea8./10^6,VolumeChange8./86400,VolumeChangeUncer8./86400,VolumeChangeUncer8./86400,SubmergedAreaUncer8./10^6,SubmergedAreaUncer8./10^6,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

xlabel('Submerged area (km^2)','fontsize',16)
ylabel('Melt flux (m^3/s)','fontsize',16)
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 7],'xtick',[0:0.5:2.5]);
legend({'20200318-20200404','20190327-20190417'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(SubmergedArea1./10^6,VolumeChange1./86400,VolumeChangeUncer1./86400,VolumeChangeUncer1./86400,SubmergedAreaUncer1./10^6,SubmergedAreaUncer1./10^6,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,2); set(gca,'fontsize',16);
errorbar(SubmergedArea3./10^6,VolumeChange3./86400,VolumeChangeUncer3./86400,VolumeChangeUncer3./86400,SubmergedAreaUncer3./10^6,SubmergedAreaUncer3./10^6,'sr','MarkerSize',5,'MarkerFaceColor','red'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 7],'xtick',[0:0.5:2.5]);
legend({'20190418-20190423','20180418-20180508'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea2./10^6,VolumeChange2./86400,VolumeChangeUncer2./86400,VolumeChangeUncer2./86400,SubmergedAreaUncer2./10^6,SubmergedAreaUncer2./10^6,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea7./10^6,VolumeChange7./86400,VolumeChangeUncer7./86400,VolumeChangeUncer7./86400,SubmergedAreaUncer7./10^6,SubmergedAreaUncer7./10^6,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 

subplot(1,4,3); set(gca,'fontsize',16);
errorbar(SubmergedArea9./10^6,VolumeChange9./86400,VolumeChangeUncer9./86400,VolumeChangeUncer9./86400,SubmergedAreaUncer9./10^6,SubmergedAreaUncer9./10^6,'sm','MarkerSize',5,'MarkerFaceColor','magenta'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 7],'xtick',[0:0.5:2.5]);
legend({'20200505-20200606','20160427-20160603', '20210525-20210605'},'Location','northwest','fontsize',14)
grid on
hold off

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(SubmergedArea5./10^6,VolumeChange5./86400,VolumeChangeUncer5./86400,VolumeChangeUncer5./86400,SubmergedAreaUncer5./10^6,SubmergedAreaUncer5./10^6,'sg','MarkerSize',5,'MarkerFaceColor','green'); hold on;

subplot(1,4,4); set(gca,'fontsize',16);
errorbar(SubmergedArea6./10^6,VolumeChange6./86400,VolumeChangeUncer6./86400,VolumeChangeUncer6./86400,SubmergedAreaUncer6./10^6,SubmergedAreaUncer6./10^6,'dg','MarkerSize',5,'MarkerFaceColor','green'); 

xlabel('Submerged area (km^2)','fontsize',16)
% ylabel('Melt flux (m^3/s)')
axis tight
set(gca,'xlim',[0 2.5],'ylim',[0 7],'xtick',[0:0.5:2.5]);
legend({'20160706-20160729','20210715-20210727'},'Location','northwest','fontsize',14)
grid on
hold off

saveas(gcf,[figure_path,'MeltFluxvsSubmergedArea.eps'],'epsc');
% sgtitle('Meltwater Flux')

%% Section 7: Testing Different plotting methods
% Following Code plots draft on Y-axis with 0 at the top and Melt rate on
% the x-axis

% figure(5), clf
% subplot(1,4,1)
% errorbar(VolumeChange4,Draft4,Draftrange4./2,Draftrange4./2,VolumeChangeUncer4,VolumeChangeUncer4,'sr','MarkerSize',5,'MarkerFaceColor','red'); hold on;
% 
% errorbar(VolumeChange8,Draft8,Draftrange8./2,Draftrange8./2,VolumeChangeUncer8,VolumeChangeUncer8,'sb','MarkerSize',5,'MarkerFaceColor','blue');
% 
% set(gca, 'YDir','reverse')
% 
% % plot(Draft2,Delt_V2,'sr','markerfacecolor','r');
% xlabel('Melt rate (m/d)')
% ylabel('Draft (m)')
% axis tight
% legend({'20200318-20200404','20190327-20190417'},'Location','northwest')
% grid on
% hold off
% 
% subplot(1,4,2)
% errorbar(VolumeChange1,Draft1,Draftrange1./2,Draftrange1./2,VolumeChangeUncer1,VolumeChangeUncer1,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;
% 
% subplot(1,4,2)
% errorbar(VolumeChange3,Draft3,Draftrange3./2,Draftrange3./2,VolumeChangeUncer3,VolumeChangeUncer3,'sr','MarkerSize',5,'MarkerFaceColor','red'); 
% 
% set(gca, 'YDir','reverse')
% 
% xlabel('Melt rate (m/d)')
% ylabel('Draft (m)')
% axis tight
% legend({'20190418-20190423','20180418-20180508'},'Location','northwest')
% grid on
% hold off
% 
% 
% subplot(1,4,3)
% errorbar(VolumeChange2,Draft2,Draftrange2./2,Draftrange2./2,VolumeChangeUncer2,VolumeChangeUncer2,'sk','MarkerSize',5,'MarkerFaceColor','black'); hold on;
% 
% subplot(1,4,3)
% errorbar(VolumeChange7,Draft7,Draftrange7./2,Draftrange7./2,VolumeChangeUncer7,VolumeChangeUncer7,'sb','MarkerSize',5,'MarkerFaceColor','blue'); 
% 
% subplot(1,4,3)
% errorbar(VolumeChange9,Draft9,Draftrange9./2,Draftrange9./2,VolumeChangeUncer9,VolumeChangeUncer9,'sm','MarkerSize',5,'MarkerFaceColor','magenta'); 
% 
% set(gca, 'YDir','reverse')
% 
% xlabel('Melt rate (m/d)')
% ylabel('Draft (m)')
% axis tight
% legend({'20200505-20200606','20160427-20160603', '20210525-20210605'},'Location','northwest')
% grid on
% hold off
% 
% 
% subplot(1,4,4)
% errorbar(VolumeChange5,Draft5,Draftrange5./2,Draftrange5./2,VolumeChangeUncer5,VolumeChangeUncer5,'sg','MarkerSize',5,'MarkerFaceColor','green'); hold on;
% 
% subplot(1,4,4)
% errorbar(VolumeChange6,Draft6,Draftrange6./2,Draftrange6./2,VolumeChangeUncer6,VolumeChangeUncer6,'dg','MarkerSize',5,'MarkerFaceColor','green'); 
% 
% set(gca, 'YDir','reverse')
% 
% xlabel('Melt rate (m/d)')
% ylabel('Draft (m)')
% axis tight
% legend({'20190418-20190423','20180418-20180508'},'Location','northwest')
% grid on
% hold off
% 


