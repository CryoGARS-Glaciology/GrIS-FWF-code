%%% create thickness maps for AeroDEM & ArcticDEM for each study site
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/';

%read in BedMachine Greenland v5 (elevations referenced to the geoid)
cd('/Users/ellynenderlin/Research/miscellaneous/BedMachine/');
bx = ncread('BedMachineGreenland-v5.nc','x');
by = ncread('BedMachineGreenland-v5.nc','y');
bz = ncread('BedMachineGreenland-v5.nc','bed'); %bed elevations
bg = ncread('BedMachineGreenland-v5.nc','geoid'); %geoid elevations relative to WGS84 ellipsoid
%figure; imagesc(bx,by,bz'); axis xy equal; colormap jet; colorbar; %plot bed to check everything looks good
%figure; imagesc(bx,by,bg'); axis xy equal; colormap jet; colorbar; %plot geoid to check everything looks good

%loop through the AeroDEMs to make site-specific thickness maps
cd([root_dir,'AeroDEM/']);
AeroDEMs = dir('*19*.tif'); %if all in the same directory
for j = 34:length(AeroDEMs)
    %load the AeroDEM
    [Z,S] = readgeoraster(AeroDEMs(j).name);
    Zo.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
    Zo.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
    Zo.z = Z; Zo.z(Zo.z==min(min(Zo.z))) = NaN; clear Z; %replace no-data value of ~-3.4028e38 with NaNs
    % figure; imagesc(Zo.x,Zo.y,Zo.z); axis xy equal; colormap jet; colorbar; %plot to check elevations look OK & no-data values are NaNs
    
    if contains(S.ProjectedCRS.Name,'UTM')
        disp(['Wrong (UTM) projection for ',AeroDEMs(j).name,' (',num2str(j),'th file)']);
        clear S Zo;
    else
    
    %create a bounding box for the DEM
    BB = [S.XWorldLimits(1),S.YWorldLimits(1); S.XWorldLimits(1),S.YWorldLimits(2);...
        S.XWorldLimits(2),S.YWorldLimits(2); S.XWorldLimits(2),S.YWorldLimits(1)];
    
    %load the ArcticDEM for the site subset
    cd('/Users/ellynenderlin/Research/miscellaneous/ArcticDEM10m/');
    info = georasterinfo('ArcticDEM_10m_Greenland_Mosaic.tif'); %pull spatial referencing info for the FULL dem
    ArcticDEMx = single([info.RasterReference.XWorldLimits(1)+0.5*info.RasterReference.CellExtentInWorldX:info.RasterReference.CellExtentInWorldX:info.RasterReference.XWorldLimits(2)-0.5*info.RasterReference.CellExtentInWorldX]);
    ArcticDEMy = single([info.RasterReference.YWorldLimits(2)-0.5*info.RasterReference.CellExtentInWorldY:-info.RasterReference.CellExtentInWorldY:info.RasterReference.YWorldLimits(1)+0.5*info.RasterReference.CellExtentInWorldY]);
    dem_xind = [find(ArcticDEMx<=min(BB(:,1)),1,'last'),find(ArcticDEMx>=max(BB(:,1)),1,'first')];
    dem_yind = [find(ArcticDEMy>=max(BB(:,2)),1,'last'),find(ArcticDEMy<=min(BB(:,2)),1,'first')];
    I = imread('ArcticDEM_10m_Greenland_Mosaic.tif','PixelRegion',{dem_yind,dem_xind}); I(I==-9999) = NaN;
    Zf.x = ArcticDEMx(min(dem_xind):max(dem_xind));
    Zf.y = ArcticDEMy(min(dem_yind):max(dem_yind));
    Zf.z = I; clear I;
    % figure; imagesc(Zf.x,Zf.y,Zf.z); axis xy equal; colormap jet; colorbar; %plot to double check
    
    %subset the bed & geoid maps and create dummy (NaN) thickness maps for each DEM
    bed_xind = [find(bx<=min(BB(:,1)),1,'last'),find(bx>=max(BB(:,1)),1,'first')];
    bed_yind = [find(by>=max(BB(:,2)),1,'last'),find(by<=min(BB(:,2)),1,'first')];
    bed_xsub = bx(min(bed_xind):max(bed_xind)); bed_ysub = by(min(bed_yind):max(bed_yind));
    bed_zsub = bz(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    bed_gsub = bg(min(bed_xind):max(bed_xind),min(bed_yind):max(bed_yind))'; %rotate so rows are y, columns are x
    % figure; imagesc(bed_xsub,bed_ysub,bed_zsub); axis xy equal; colormap jet; colorbar; %plot to check subsetting
    Ho = NaN(size(bed_zsub)); Hf = NaN(size(bed_zsub));
    
    %interpolate the surface elevations to the subset geoid map
    [zo_xgrid,zo_ygrid] = meshgrid(Zo.x,Zo.y);
    [zf_xgrid,zf_ygrid] = meshgrid(Zf.x,Zf.y);
    [b_xgrid,b_ygrid] = meshgrid(bed_xsub,bed_ysub);
    zo_interp = single(interp2(zo_xgrid,zo_ygrid,double(Zo.z),single(b_xgrid),single(b_ygrid),'cubic'));
    zf_interp = single(interp2(zf_xgrid,zf_ygrid,double(Zf.z),single(b_xgrid),single(b_ygrid),'cubic'));
    
    %convert interpolated ellipsoidal elevations to orthometric elevations
    zo_interp = zo_interp - single(bed_gsub); zf_interp = zf_interp - single(bed_gsub);
    
    %calculate the ice thickness for each DEM under the assumption of flotation
    floating = 1026/(1026-900); %multiplier to convert elevation to thickness
    Ho_float = floating.*zo_interp; Hf_float = floating.*zf_interp;
    
    %calculate the ice thickness within each DEM as the difference between
    %the surface & bed elevations
    Ho_ground = zo_interp - bed_zsub; Hf_ground = zf_interp - bed_zsub;

    %if the floating thickness is less than the difference thickness, fill
    %in the dummy thickness map with the floating thickness, else use the
    %difference thickness
    Ho(Ho_float<Ho_ground) = Ho_float(Ho_float<Ho_ground);
    Ho(Ho_float>=Ho_ground) = Ho_ground(Ho_float>=Ho_ground);
    Hf(Hf_float<Hf_ground) = Hf_float(Hf_float<Hf_ground);
    Hf(Hf_float>=Hf_ground) = Hf_ground(Hf_float>=Hf_ground);
    Ho(Ho<0) = 0; Hf(Hf<0) = 0;
    %check the thickness maps
    figure; set(gca,'position',[50 50 800 400]); 
    cmap = colormap(jet(10001)); cmap(1,:) = [1 1 1];
    clims = [-1, max([max(Ho(~isnan(Ho))),max(Hf(~isnan(Hf)))])];
    subplot(1,2,1);
    imagesc(bed_xsub,bed_ysub,Ho); axis xy equal; colormap(gca,cmap);
    set(gca,'clim',clims,'fontsize',12);
    xlabel('Easting (m)','fontsize',12); ylabel('Northing (m)','fontsize',12); 
    subplot(1,2,2); pos = get(gca,'position');
    imagesc(bed_xsub,bed_ysub,Hf); axis xy equal; colormap(gca,cmap); cbar = colorbar;
    set(gca,'clim',clims,'yticklabel',[],'fontsize',12); cbar.Label.String = 'thickness (m)';
    xlabel('Easting (m)','fontsize',12); 
    set(gca,'position',[pos(1)-0.05 pos(2) pos(3) pos(4)]);
    saveas(gcf,[root_dir,'terminus-thickness/',AeroDEMs(j).name(1:end-8),'thickness_maps.png'],'png');
    
    %save the maps as geotiffs
    cd([root_dir,'terminus-thickness/']);
    So = S;
    So.XWorldLimits = sort([double(nanmean(bx(min(bed_xind)-1:min(bed_xind)))) double(nanmean(bx(max(bed_xind):max(bed_xind)+1)))]);
    So.YWorldLimits = sort([double(nanmean(by(min(bed_yind)-1:min(bed_yind)))) double(nanmean(by(max(bed_yind):max(bed_yind)+1)))]);
    So.RasterSize = size(Ho); 
    geotiffwrite([AeroDEMs(j).name(1:end-8),'thickness_',AeroDEMs(j).name(end-7:end-4),'.tif'],Ho,So,'CoordRefSysCode',3413);
    Sf = S;
    Sf.XWorldLimits = sort([double(nanmean(bx(min(bed_xind)-1:min(bed_xind)))) double(nanmean(bx(max(bed_xind):max(bed_xind)+1)))]);
    Sf.YWorldLimits = sort([double(nanmean(by(min(bed_yind)-1:min(bed_yind)))) double(nanmean(by(max(bed_yind):max(bed_yind)+1)))]);
    Sf.RasterSize = size(Hf); 
    geotiffwrite([AeroDEMs(j).name(1:end-8),'thickness_2016.tif'],Hf,Sf,'CoordRefSysCode',3413);
    
    %clear out variables
    clear Z* S* BB ArcticDEM* dem*ind bed_* b_*grid H* zo* zf* clims;
    end
    cd([root_dir,'AeroDEM/']);
end