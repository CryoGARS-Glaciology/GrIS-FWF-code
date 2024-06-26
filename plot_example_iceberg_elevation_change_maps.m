%%% Plot example iceberg elevation change figures
clearvars; close all;

addpath('/users/ellynenderlin/Research/miscellaneous/general-code','/users/ellynenderlin/Research/miscellaneous/general-code/cmocean');

% data_path = '/Users/ellynenderlin/Research/NSF_Greenland-Calving/iceberg-calving/JI/';
data_path = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-meltrates/';
figure_path = ['/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/figures/iceberg-melt/'];
cd(data_path);

%specify plot params
cmax = 70; %saturated elevation in m a.s.l.
zoom_buff = 1000; %+/- m for buffer around iceberg
colors = 'thermal'; %cmocean color pallette
fonts = 24;

%load and make map figures of each DEM
% DEMs = dir('JI*melange-DEM.mat');
DEMs = dir('H*DEM.mat');
for j = 1:2 %1:length(DEMs)
    %load the data
    load(DEMs(j).name);
    if exist('Y') == 1; Z = Y; clear Y; end
    
    %replace negative elevations with zeros for visualization
    if isfield(Z,'z_elpsd_adjust')
        z = Z.z_elpsd_adjust;
    else
        z = Z.z.ortho;
    end
    z(z<0) = 0;
    
    %plot
    figure; set(gcf,'position',[50 50 1200 800]);
    imagesc(Z.x,Z.y,z-nanmedian(min(z))); axis xy equal; hold on;
    set(gca,'clim',[-1 cmax]);
    cmap = cmocean(colors,1001); cmap(1,:) = [1 1 1]; colormap(gca,cmap);
    cbar = colorbar; cbar.Label.String = 'elevation (m a.s.l.)';
    
    %zoom to focus on the region near the terminus
%     set(gca,'xlim',[-1.94e5 -1.84e5],'ylim',[-2.277e6 -2.269e6],...
%         'xtick',[-1.94e5:0.01e5:-1.84e5],'xticklabel',[-194:1:-184],...
%         'ytick',[-2.277e6:0.01e5:-2.269e6],'yticklabel',[-2277:1:-2269]);
    set(gca,'xlim',[3.09e5 3.21e5],'ylim',[-2.582e6 -2.575e6],...
        'xtick',[-3.09e5:0.01e5:3.21e5],'xticklabel',[309:1:321],...
        'ytick',[-2.582e6:0.01e5:-2.575e6],'yticklabel',[-2582:1:-2575]);
    set(gca,'fontsize',fonts);
    xlabel('Easting (km)','fontsize',fonts); ylabel('Northing (km)','fontsize',fonts);
%     title([DEMs(j).name(4:7),'/',DEMs(j).name(8:9),'/',DEMs(j).name(10:11)]);
    title([DEMs(j).name(5:8),'/',DEMs(j).name(9:10),'/',DEMs(j).name(11:12)]);
    drawnow;
    
    %prompt if you want to make a cropped, zoomed image with contours
    answer = questdlg('Create a zoomed map?','Zoom map?','Yes','No','No');
    switch answer
        case 'Yes'
            %identify zoom location
            disp('click on a point in the middle of an iceberg');
            [b] = ginput(1);
            
            %crop DEM
            xcrops = [find(Z.x<=b(1)-zoom_buff,1,'last'),find(Z.x>=b(1)+zoom_buff,1,'first')];
            ycrops = [find(Z.y<=b(2)+zoom_buff,1,'first'),find(Z.y>=b(2)-zoom_buff,1,'last')];
            if max(ycrops) == length(Z.y)
                ycrops = [find(Z.y>=b(2)-zoom_buff,1,'first'),find(Z.y<=b(2)+zoom_buff,1,'last')];
            end
            xcropped = Z.x(xcrops(1):xcrops(2));
            ycropped = Z.y(ycrops(1):ycrops(2));
            zcropped = z(ycrops(1):ycrops(2),xcrops(1):xcrops(2));
            
            %plot cropped figure
            figure; set(gcf,'position',[50 50 1200 800]);
            imagesc(xcropped,ycropped,zcropped-nanmedian(min(z))); axis xy equal; hold on;
            set(gca,'clim',[-1 cmax]);
            cmap = cmocean(colors,1001); cmap(1,:) = [1 1 1]; colormap(gca,cmap);
            cbar = colorbar; cbar.Label.String = 'elevation (m a.s.l.)';
            
            %format
            set(gca,'xlim',[min(xcropped),max(xcropped)],'ylim',[min(ycropped),max(ycropped)],...
                'xtick',[ceil(min(xcropped)/100)*100:400:floor(max(xcropped)/100)*100],...
                'xticklabel',[ceil(min(xcropped)/100)/10:0.4:floor(max(xcropped)/100)/10],...
                'ytick',[ceil(min(ycropped)/100)*100:400:floor(max(ycropped)/100)*100],...
                'yticklabel',[ceil(min(ycropped)/100)/10:0.4:floor(max(ycropped)/100)/10]);
            set(gca,'fontsize',fonts);
            xlabel('Easting (km)','fontsize',fonts); ylabel('Northing (km)','fontsize',fonts);
%             title([DEMs(j).name(4:7),'/',DEMs(j).name(8:9),'/',DEMs(j).name(10:11)]);
            title([DEMs(j).name(5:8),'/',DEMs(j).name(9:10),'/',DEMs(j).name(11:12)]);
            
            %add contour lines
            [cont,conth] = contour(xcropped,ycropped,zcropped,[0:2:cmax+20]);
            conth.LineColor = 'k';
            drawnow;
            
            %save zoomed map
            saveas(gcf,[figure_path,DEMs(j).name(1:11),'_zoomed-iceberg-elevation-map.eps'],'epsc');
            saveas(gcf,[figure_path,DEMs(j).name(1:11),'_zoomed-iceberg-elevation-map.png'],'png');
            
            clear *crops *cropped cont*;
        case 'No'
            disp('moving on...');
    end
    

    
    %clear data
    clear Z z;
end