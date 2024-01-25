%PLOT_MELTRATE_FIGS_UPDATED: Generate maps and melt rate subplots for all
%Antarctic iceberg melt datasets. Includes ocean observation data in
%figures. Exports concatenated data for all sites to a table.


%% Initialize
clearvars; close all; drawnow;
addpath('/users/ellynenderlin/Research/miscellaneous/general-code','/users/ellynenderlin/Research/miscellaneous/general-code/cmocean');
% 
% %specify paths & file names for data
% iceberg_path = '/Users/adamdickson/Desktop/CSVFilesIcebergMelt';
% figure_path = [iceberg_path,'figures/'];
iceberg_path = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-meltrates/';
figure_path = ['/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/figures/iceberg-melt/'];

%specify generic variables
plot_yrs = []; avgx = []; avgy = []; region = []; warning off;
rho_sw = 1026; %sea water density in kg m^-3
years = [2011.75 2022.25]; year_ticks = [2013:2:2022]; %approximate date range for plots
plot_marker = 's';
symbol_size = 8;


%% create subplots & add site info to the map
close all;

%set-up subplots for graphs
figureB = figure; set(gcf,'position',[50 400 400 400]); %all data lumped
% sub1b = subplot(1,2,1); sub2b = subplot(1,2,2); %meltwater flux on left, meltrate on right
figureC = figure; set(gcf,'position',[450 400 800 1200]); %subplot for each region

cd(iceberg_path);
disp('Generating plots...');

%set-up dummy vectors to fill with concatenated variables for all sites
start_yr = []; end_yr = []; avg_x = []; avg_y = []; depth = []; depth_uncert = []; subarea = []; subarea_uncert = [];
meltflux = []; meltflux_uncert = []; meltrate = []; meltrate_uncert = [];
%start & end coordinates for each iceberg
xcoord_o = []; ycoord_o = [];
xcoord_f = []; ycoord_f = [];
region_flag = [];

% Finds all meltinfo.csv files 
MeltFiles = dir('*meltinfo.csv');
% Grabs site abbreviation for each meltinfo.csv file
for g = 1:length(MeltFiles)
    underlineLocation = strfind(MeltFiles(g).name, '_');
    site(g,:) = MeltFiles(g).name(1:underlineLocation-1);
end

%load site and region list
glacial_abbrev_region = readtable('GlacialAbbreWithRegionlocation');
abbrev = glacial_abbrev_region(:,1); abbrev = table2array(abbrev);
glacial_regions = glacial_abbrev_region{:,2};
abbrevs = string(abbrev);
clear abbrev glacial_abbrev_region;

%identify unique regions, assigns indices & colors
% regions = string(unique(glacial_regions));
regions = ["NW";"NO";"CW";"CE";"SW";"SE"]; %NEED TO CHANGE NO TO NE IF USING NE FOR ZIM + SUQ IN REGION TABLE
% cmcolor = cmocean('phase',length(regions)+1); cmcolor = cmcolor(1:end-1,:);
cmcolor = [118,42,131; 27,120,55; 153,112,171; 90,174,97; 194,165,207; 166,219,160]./255;

%create dummy matrices to concatenate data
avgx = []; avgy = []; meltrate_v_draft = []; 
%     date_o = []; xcoord_o = []; ycoord_o = [];
%     date_f = []; xcoord_f = []; ycoord_f = [];
%     flux = []; sub_area = []; meltrate = []; keeld = [];
melt = struct('draft',[],'Asub',[],'m',[],'dVdt',[]);
for j = 2:length(regions); melt(j) = melt(1); end

%loop through each date pair & plot data
for j = 1:length(MeltFiles)
    M=readtable(MeltFiles(j).name); %table exported using plot_export_iceberg_melt_data.m
    %         disp(M.Properties.VariableNames); %uncomment to display table headers prior to array conversion
    M = table2array(M); %convert to array to enable easier indexing
    disp(['date: ',MeltFiles(j).name(5:12),'-',MeltFiles(j).name(14:21)]);
    site = MeltFiles(j).name(1:3); disp(site);

    %find the index for the site name in the list of names and regions
    for k = 1:length(abbrevs)
        if strcmp(abbrevs(k,:),site) == 1
            site_ref = k;
        end
    end

    %get the region index (and therefore color) for the site
    for k = 1:length(regions)
        if strcmp(regions(k,:),string(glacial_regions(site_ref,:))) == 1
            region_ref = k;
        end
    end
    disp(regions(region_ref));
    
    
    %identify data with clear issues
    bad_data = find(M(:,18)<0); M(bad_data,:) = [];
    
    %pull variables
    dt = M(:,1); %time
    xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); Vo = M(:,6); %initial locations, median elev, density, volume
    xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); Vf = M(:,11); %same as above but final
    coregzo = M(:,12); coregzf = M(:,13);
    dz = M(:,14); dz_sigma = M(:,15);
    dVdt = M(:,16); dVdt_uncert = M(:,17);
    
    %recalculate draft & submerged area to make sure they are consistent (methods may have been adjusted slightly over time)
    draft = (nanmean([rhoo rhof],2)./(repmat(rho_sw,length(nanmean([rhoo rhof],2)),1)-nanmean([rhoo rhof],2))).*nanmean([zo zf],2); %draft = M(:,18);
    draft_uncert = M(:,19);
    Asurf = M(:,20); Asurf_uncert = M(:,21);
    lat_area = M(:,22) - Asurf; perim = lat_area./draft; clear lat_area; lat_area = perim.*draft; Asub = lat_area + Asurf; clear lat_area perim; %Asub = M(:,22);
    Asub_uncert = M(:,23);
    m = dVdt./Asub; %melt rate variable for plotting
    disp(['average increase in melt rate with draft: ',num2str(round(nanmean(365*m./draft),4)),' m/yr per m depth']);
    
    %         date_o = [date_o; repmat(str2num(meltinfo(j).name(4:11)),size(xo))]; xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo];
    %         date_f = [date_f; repmat(str2num(meltinfo(j).name(13:20)),size(xf))]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf];
    %         flux = [flux; dVdt]; sub_area = [sub_area; Asub]; meltrate = [meltrate; (dVdt./Asub)]; keeld = [keeld; draft];
    
    %compile data to create a concatenated table for all sites, all dates
    decidate_o = convert_to_decimaldate(MeltFiles(j).name(5:12));
    decidate_f = convert_to_decimaldate(MeltFiles(j).name(14:21));
    start_yr = [start_yr; repmat(decidate_o,length(draft),1)]; end_yr = [end_yr; repmat(decidate_f,length(draft),1)];
    plot_yrs = [decidate_o decidate_f];
    clear decidate_*;
    avgx = [avgx; nanmean([xo xf],2)]; avgy = [avgy; nanmean([yo yf],2)]; %average coordinates for regional map
    xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf]; %coordinates for coordinate data table
    avg_x = [avg_x; nanmean([xo xf],2)]; avg_y = [avg_y; nanmean([yo yf],2)]; %average coordinates for site map & regional data table
    depth = [depth; draft]; depth_uncert = [depth_uncert; draft_uncert]; %keel depth
    subarea = [subarea; Asub]; subarea_uncert = [subarea_uncert; Asub_uncert]; %submerged area
    meltflux = [meltflux; dVdt]; meltflux_uncert = [meltflux_uncert; dVdt_uncert]; %meltwater flux
    meltrate = [meltrate; m]; meltrate_uncert = [meltrate_uncert; abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2)]; %melt rate
    
    %compile region refs
    region_flag = [region_flag; repmat(region_ref,size(m))];
    
    %add to a structure for each region
    melt(region_ref).draft = [melt(region_ref).draft; draft];
    melt(region_ref).Asub = [melt(region_ref).Asub; Asub];
    melt(region_ref).dVdt = [melt(region_ref).dVdt; dVdt];
    melt(region_ref).m = [melt(region_ref).m; m];
    
    %display size range data
    disp(['Surface area range: ',num2str(min(Asurf)),' - ',num2str(max(Asurf)),' m^2']);
    disp(['Draft range: ',num2str(min(draft)),' - ',num2str(max(draft)),' m']);
    disp(['Submerged area range: ',num2str(min(Asub)),' - ',num2str(max(Asub)),' m^2']);
    
    %multi-panel subplots of all data
    figure(figureB);
%     subplot(sub1b);
    errorbar(Asub,dVdt/86400,dVdt_uncert/86400,dVdt_uncert/86400,Asub_uncert,Asub_uncert,plot_marker,...
        'color','k','markerfacecolor',cmcolor(region_ref,:),'markersize',symbol_size./2,'markeredgecolor',cmcolor(region_ref,:)); hold on;
    plot(Asub,dVdt/86400,plot_marker,'markerfacecolor',cmcolor(region_ref,:),'markeredgecolor',cmcolor(region_ref,:),...
        'markersize',symbol_size./2,'linewidth',1); hold on;
%     subplot(sub2b);
%     errorbar(draft,m,abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
%         'color','k','markerfacecolor',cmcolor(region_ref,:),'markersize',symbol_size,'markeredgecolor','k'); hold on;
        
    % [f,gof] = fit(draft,m,'poly1'); %fit a linear trendline
    % ci = confint(f,0.95); % 95% confidence interval
    % % plot the line
    % plot([0:10:max(draft)],...
    %     feval(f,[0:10:max(draft)]),...
    %     '-','color','r','linewidth',2); hold on;
    % % plot the confidence interval as a semi-transparent polygon
    % fill([0:10:max(draft),max(draft):-10:0],...
    %     [(ci(1,1).*[0:10:max(draft)]+ci(1,2)),...
    %     (ci(2,1).*[max(draft):-10:0]+ci(2,2))],...
    % 'r','FaceAlpha',0.25,'EdgeColor','r');

    
    %regional subplots
    figure(figureC);
    subplot(ceil(length(regions)/2),2,region_ref); %2 columns for west vs east coast
    errorbar(draft,m,abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
        'color','k','markerfacecolor',cmcolor(region_ref,:),'markersize',symbol_size./2,'markeredgecolor',cmcolor(region_ref,:)); hold on;
    
    % 
    %remove date-specific variables
    clear site region_ref site_ref;
    % clear m dt xo yo zo rhoo Vo xf yf zf rhof Vf coregzo coregzf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
end
%format lumped figure
figure(figureB);
% subplot(sub1b); 
grid on; ylims = get(gca,'ylim'); xlims = get(gca,'xlim'); xticks = get(gca,'xtick');
set(gca,'ylim',[0 18],'xlim',[0 max(xlims)],'xtick',xticks,'xticklabel',xticks/10^6,'fontsize',20);
xlabel('Submerged area (km^2)','fontsize',20); ylabel('Meltwater flux (m^3/s)','fontsize',20);
% subplot(sub2b); 
% grid on; ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
% set(gca,'ylim',[0 1],'xlim',[0 450],'fontsize',20);
% xlabel('Draft (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);


%formt regional subplots
figure(figureC);
for j = 1:length(regions)
    subplot(ceil(length(regions)/2),2,j); %2 columns for west vs east coast
    grid on; 
    set(gca,'ylim',[0 1.0],'xlim',[0 500],'fontsize',16);
    ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
%     title(regions(j));
    
    
    %adjust subplot positions
%     pos = get(gca,'position');
% %     set(gca,'position',[pos(1) pos(2)+0.02 pos(3) pos(4)]);
%     set(gca,'position',[pos(1) pos(2) 1.05*pos(3) 1.05*pos(4)]);
    

end



%formt regional subplots
figure(figureC);
for j = 1:length(regions)
    subplot(ceil(length(regions)/2),2,j);
    grid on; ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
    set(gca,'ylim',[0 1],'xlim',[0 500],'ytick',[0:0.2:1],'fontsize',16);
    text(0.075*max(xlims),0.90*max(ylims),regions(j),'fontsize',20);
    
    %add trendlines
    [f,gof] = fit(melt(j).draft,melt(j).m,'poly1'); %fit a linear trendline
    ci = confint(f,0.95); % 95% confidence interval
    % plot the line
    plot([0:10:max(melt(j).draft)],...
        feval(f,[0:10:max(melt(j).draft)]),...
        '-','color',cmcolor(j,:),'linewidth',2); hold on;
    % plot the confidence interval as a semi-transparent polygon
    fill([0:10:max(melt(j).draft),max(melt(j).draft):-10:0],...
        [(ci(1,1).*[0:10:max(melt(j).draft)]+ci(1,2)),...
        (ci(2,1).*[max(melt(j).draft):-10:0]+ci(2,2))],...
    cmcolor(j,:),'FaceAlpha',0.25,'EdgeColor',cmcolor(j,:));


    %add x axis labels along the bottom row
    if ceil(j/2) == ceil(length(regions)/2)
        xlabel('Draft (m)','fontsize',16); 
    end
    
    %add y axis labels along the left column
    if mod(j,2) == 1
        ylabel('Melt rate (m/d)','fontsize',16);
    end
end
% M(j).draft, M(j).m
% Trying to add tend line
% [f,gof] = fit(draft,m,'poly1'); %fit a linear trendline
% ci = confint(f,0.95); % 95% confidence interval
% % plot the line
% plot([0:10:max(draft)],...
%     feval(f,[0:10:max(draft)]),...
%     '-','color','r','linewidth',2); hold on;
% % plot the confidence interval as a semi-transparent polygon
% fill([0:10:max(draft),max(draft):-10:0],...
%     [(ci(1,1).*[0:10:max(draft)]+ci(1,2)),...
%     (ci(2,1).*[max(draft):-10:0]+ci(2,2))],...
% 'r','FaceAlpha',0.25,'EdgeColor','r');
    
%add axis labels to the bottom plot
% AddLetters2Plots(figureC, {'(a)', '(b)', '(c)' , '(d)', '(e)', '(f)'}, 'Direction', 'TopDown')
% GLI=axes(figureC,'visible','off'); 
% GLI.Title.Visible='on';
% GLI.XLabel.Visible='on';
% GLI.YLabel.Visible='on';
% xlabel(GLI,'Draft (m)','fontsize',20);
% ylabel(GLI,{'Melt rate (m/d)';''},'fontsize',20);
% title(GLI,{'Regional Melt Rate (m/d)';''},'fontsize',20);

%save the figures
saveas(figureB,[figure_path,'GrIS_regional-meltflux-vs-area_plot.eps'],'epsc');
saveas(figureC,[figure_path,'GrIS_regional-meltrate-vs-draft_subplots.eps'],'epsc');


%% regional subplots with months shown in colors
close all;

%set-up subplots for graphs
figureB = figure; set(gcf,'position',[50 400 800 400]); %all data lumped
sub1b = subplot(1,2,1); sub2b = subplot(1,2,2); %meltwater flux on left, meltrate on right
figureC = figure; set(gcf,'position',[450 400 800 1200]); %subplot for each region
% mo_cmap = [cmocean('thermal',6); flipud(cmocean('haline',6))];
mo_cmap = flipud([80,0,37; 165,0,38; 215,48,39; 244,109,67; 253,174,97; 254,224,144;...
224,243,248; 171,217,233; 116,173,209; 69,117,180; 49,54,149; 30,20,100])./255;

% cd('/Users/adamdickson/Desktop/CSVFilesIcebergMelt');
cd('/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-meltrates/');
disp('Generating plots...');

%set-up dummy vectors to fill with concatenated variables for all sites
start_yr = []; end_yr = []; avg_x = []; avg_y = []; depth = []; depth_uncert = []; subarea = []; subarea_uncert = [];
meltflux = []; meltflux_uncert = []; meltrate = []; meltrate_uncert = [];
%start & end coordinates for each iceberg
xcoord_o = []; ycoord_o = [];
xcoord_f = []; ycoord_f = [];

% Finds all meltinfo.csv files 
MeltFiles = dir('*meltinfo.csv');
% Grabs site abbreviation for each meltinfo.csv file
for g = 1:length(MeltFiles)
    underlineLocation = strfind(MeltFiles(g).name, '_');
    site(g,:) = MeltFiles(g).name(1:underlineLocation(1)-1);
    decidate(1) = convert_to_decimaldate(MeltFiles(g).name(underlineLocation(1)+1:underlineLocation(1)+8));
    decidate(2) = convert_to_decimaldate(MeltFiles(g).name(underlineLocation(1)+10:underlineLocation(1)+17));
    sitedate(g,:) = nanmean([decidate(1),decidate(2)]);
    sitemo(g,:) = round(12*nanmean([(decidate(1)-floor(decidate(1))),(decidate(2)-floor(decidate(2)))]));
end

%load site and region list
glacial_abbrev_region = readtable('GlacialAbbreWithRegionlocation');
abbrev = glacial_abbrev_region(:,1); abbrev = table2array(abbrev);
glacial_regions = glacial_abbrev_region{:,2};
abbrevs = string(abbrev);
clear abbrev glacial_abbrev_region;

%identify unique regions, assigns indices & colors
% regions = string(unique(glacial_regions));
regions = ["NW";"NO";"CW";"CE";"SW";"SE"]; %NEED TO CHANGE NO TO NE IF USING NE FOR ZIM + SUQ IN REGION TABLE

%plot
avgx = []; avgy = []; meltrate_v_draft = [];
%     date_o = []; xcoord_o = []; ycoord_o = [];
%     date_f = []; xcoord_f = []; ycoord_f = [];
%     flux = []; sub_area = []; meltrate = []; keeld = [];

%loop through each date pair & plot data
for j = 1:length(MeltFiles)
    M=readtable(MeltFiles(j).name); %table exported using plot_export_iceberg_melt_data.m
    %         disp(M.Properties.VariableNames); %uncomment to display table headers prior to array conversion
    M = table2array(M); %convert to array to enable easier indexing
    disp(['date: ',MeltFiles(j).name(5:12),'-',MeltFiles(j).name(14:21)]);
    site = MeltFiles(j).name(1:3); disp(site);
    
    %find the index for the site name in the list of names and regions
    for k = 1:length(abbrevs)
        if strcmp(abbrevs(k,:),site) == 1
            site_ref = k;
        end
    end

    %get the region index (and therefore color) for the site
    for k = 1:length(regions)
        if strcmp(regions(k,:),string(glacial_regions(site_ref,:))) == 1
            region_ref = k;
        end
    end
    disp(regions(region_ref));

    %identify data with clear issues
    bad_data = find(M(:,18)<0); M(bad_data,:) = [];

    %pull variables
    dt = M(:,1); %time
    xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); Vo = M(:,6); %initial locations, median elev, density, volume
    xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); Vf = M(:,11); %same as above but final
    coregzo = M(:,12); coregzf = M(:,13);
    dz = M(:,14); dz_sigma = M(:,15);
    dVdt = M(:,16); dVdt_uncert = M(:,17);
    
    %recalculate draft & submerged area to make sure they are consistent (methods may have been adjusted slightly over time)
    draft = (nanmean([rhoo rhof],2)./(repmat(rho_sw,length(nanmean([rhoo rhof],2)),1)-nanmean([rhoo rhof],2))).*nanmean([zo zf],2); %draft = M(:,18);
    draft_uncert = M(:,19);
    Asurf = M(:,20); Asurf_uncert = M(:,21);
    lat_area = M(:,22) - Asurf; perim = lat_area./draft; clear lat_area; lat_area = perim.*draft; Asub = lat_area + Asurf; clear lat_area perim; %Asub = M(:,22);
    Asub_uncert = M(:,23);
    m = dVdt./Asub; %melt rate variable for plotting
    disp(['average increase in melt rate with draft: ',num2str(round(nanmean(365*m./draft),4)),' m/yr per m depth']);
    
    %         date_o = [date_o; repmat(str2num(meltinfo(j).name(4:11)),size(xo))]; xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo];
    %         date_f = [date_f; repmat(str2num(meltinfo(j).name(13:20)),size(xf))]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf];
    %         flux = [flux; dVdt]; sub_area = [sub_area; Asub]; meltrate = [meltrate; (dVdt./Asub)]; keeld = [keeld; draft];
    
    %compile data to create a concatenated table for all sites, all dates
    decidate_o = convert_to_decimaldate(MeltFiles(j).name(5:12));
    decidate_f = convert_to_decimaldate(MeltFiles(j).name(14:21));
    start_yr = [start_yr; repmat(decidate_o,length(draft),1)]; end_yr = [end_yr; repmat(decidate_f,length(draft),1)];
    plot_yrs = [decidate_o decidate_f];
    clear decidate_*;
    avgx = [avgx; nanmean([xo xf],2)]; avgy = [avgy; nanmean([yo yf],2)]; %average coordinates for regional map
    xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf]; %coordinates for coordinate data table
    avg_x = [avg_x; nanmean([xo xf],2)]; avg_y = [avg_y; nanmean([yo yf],2)]; %average coordinates for site map & regional data table
    depth = [depth; draft]; depth_uncert = [depth_uncert; draft_uncert]; %keel depth
    subarea = [subarea; Asub]; subarea_uncert = [subarea_uncert; Asub_uncert]; %submerged area
    meltflux = [meltflux; dVdt]; meltflux_uncert = [meltflux_uncert; dVdt_uncert]; %meltwater flux
    meltrate = [meltrate; m]; meltrate_uncert = [meltrate_uncert; abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2)]; %melt rate
    
    %display size range data
    disp(['Surface area range: ',num2str(min(Asurf)),' - ',num2str(max(Asurf)),' m^2']);
    disp(['Draft range: ',num2str(min(draft)),' - ',num2str(max(draft)),' m']);
    disp(['Submerged area range: ',num2str(min(Asub)),' - ',num2str(max(Asub)),' m^2']);
    
    %multi-panel subplots of all data
    figure(figureB);
    subplot(sub1b);
%     errorbar(Asub,dVdt/86400,dVdt_uncert/86400,dVdt_uncert/86400,Asub_uncert,Asub_uncert,plot_marker,...
%         'color','k','markerfacecolor',mo_cmap(sitemo(j,:),:),'markersize',symbol_size,'markeredgecolor','k'); hold on;
    plot(Asub,dVdt/86400,'+','markerfacecolor',mo_cmap(sitemo(j,:),:),...
        'markersize',symbol_size./2,'linewidth',1); hold on;
    subplot(sub2b);
%     errorbar(draft,m,abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
%         'color','k','markerfacecolor',mo_cmap(sitemo(j,:),:),'markersize',symbol_size,'markeredgecolor','k'); hold on;
    plot(draft,m,'+','markerfacecolor',mo_cmap(sitemo(j,:),:),...
        'markersize',symbol_size./2,'linewidth',1); hold on;

    %regional subplots
    figure(figureC);
    subplot(ceil(length(regions)/2),2,region_ref); %2 columns for west vs east coast
    errorbar(draft,m,abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
        'color',mo_cmap(sitemo(j,:),:),'markerfacecolor',mo_cmap(sitemo(j,:),:),'markersize',symbol_size./2,'markeredgecolor','k'); hold on;
   
    %remove date-specific variables
    clear site region_ref site_ref;
    clear m dt xo yo zo rhoo Vo xf yf zf rhof Vf coregzo coregzf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
end
%format lumped figure
figure(figureB);
subplot(sub1b); grid on; ylims = get(gca,'ylim'); xlims = get(gca,'xlim'); xticks = get(gca,'xtick');
set(gca,'ylim',[0 max(ylims)],'xlim',[0 max(xlims)],'xtick',xticks,'xticklabel',xticks/10^6,'fontsize',16);
xlabel('Submerged area (km^2)','fontsize',16); ylabel('Meltwater flux (m^3/s)','fontsize',16);
subplot(sub2b); grid on; ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
set(gca,'ylim',[0 max(ylims)],'xlim',[0 max(xlims)],'fontsize',16);
xlabel('Draft (m)','fontsize',16); ylabel('Melt rate (m/d)','fontsize',16);

%formt regional subplots
figure(figureC);
for j = 1:length(regions)
    subplot(ceil(length(regions)/2),2,j); %2 columns for west vs east coast
    grid on; 
    set(gca,'ylim',[0 1.0],'xlim',[0 500],'fontsize',16);
    ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
%     title(regions(j));
    text(0.075*max(xlims),0.90*max(ylims),regions(j),'fontsize',20);
    
    %adjust subplot positions
%     pos = get(gca,'position');
% %     set(gca,'position',[pos(1) pos(2)+0.02 pos(3) pos(4)]);
%     set(gca,'position',[pos(1) pos(2) 1.05*pos(3) 1.05*pos(4)]);
    
    %add x axis labels along the bottom row
    if ceil(j/2) == ceil(length(regions)/2)
        xlabel('Draft (m)','fontsize',16); 
    end
    
    %add y axis labels along the left column
    if mod(j,2) == 1
        ylabel('Melt rate (m/d)','fontsize',16);
    end
end
subplot(ceil(length(regions)/2),2,ceil(length(regions)/2));
for k = 1:length(mo_cmap)
    annotation('line',[0.2+0.05*k 0.2+0.05*k],[0.03 0.06],'color',mo_cmap(k,:),'linewidth',40);
end
annotation('textbox',[0.225 0.01 0.7 0.05],...
    'string','Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec','fontsize',16,...
    'edgecolor','none');

%save the figures
saveas(figureB,[figure_path,'GrIS_meltflux_meltrate_lumped_scatterplots.eps'],'epsc');
saveas(figureC,[figure_path,'GrIS_meltflux_meltrate_regional_scatterplots.eps'],'epsc');


%% save the data tables (MODIFY FOR GREENLAND)
%MELT
column_names = {'Start Date' 'End Date' 'Polar Stereo Easting' 'Polar Stereo Northing'...
    'Average Median Draft' 'Median Draft Variability' 'Average Submerged Area' 'Submerged Area Variability'...
    'Meltwater Flux' 'Meltwater Flux Uncertainty' 'Melt Rate' 'Melt Rate Uncertainty'};
column_units = {'years' 'years' 'meters' 'meters'...
    'meters b.s.l.' 'meters b.s.l.' 'cubic meters' 'cubic meters'...
    'cubic meters per day' 'cubic meters per day' 'meters per day' 'meters per day'};
T=table(start_yr,end_yr,avg_x,avg_y,depth,depth_uncert,subarea,subarea_uncert,meltflux,meltflux_uncert,meltrate,meltrate_uncert);
T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
writetable(T,[iceberg_path,'GrIS-iceberg-meltdata.csv']);
disp('Greenland iceberg melt rate text file written');
clear column_*;
%COORDINATES ONLY
column_names = {'Start Date' 'Start Polar Stereo Easting' 'Start Polar Stereo Northing'...
    'End Date' 'End Polar Stereo Easting' 'End Polar Stereo Northing' 'Region'};
column_units = {'years' 'meters' 'meters'...
    'years' 'meters' 'meters' 'unitless'};
t=table(start_yr,xcoord_o,ycoord_o,end_yr,xcoord_f,ycoord_f,region_flag);
t.Properties.VariableNames = column_names; t.Properties.VariableUnits = column_units;
writetable(t,[iceberg_path,'GrIS-iceberg-PScoords.csv']);
disp('Greenland iceberg coordinates text file written');
clear T t;


