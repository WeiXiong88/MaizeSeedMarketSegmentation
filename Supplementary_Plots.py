from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import geopandas as gpd
import rasterio
from sklearn.metrics import r2_score
from statistics import mode
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize, LogNorm
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch
#%matplotlib inline

##########################################################FOR SUPPLEMENTARY FIGURE S1#######################################
def Supplementary_Plot_S1():
    in_dir="D:\\works\\AfricaMzSg\\input\\"
    #Define color
    #type_colors = [{-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"},{-1: "red", 0: "lightgrey", 1: "green"}] #blue, red, green, maturity type
    #change_colors = {-1: "red", 0: "lightgrey", 1: "blue"} #change
    #African country boundary
    africa=gpd.read_file("E:/RSG Dropbox/Wei Xiong/Works/CurrentProcessing/0_AfricanMaizeSorghum/Plots/CountryBoundaryLines/ne_110m_admin_0_countries.shp")
    africa=africa[africa['CONTINENT']=='Africa']
    tfile=["sow_month","har_area"]
    tif_dir="D:\\works\\AfricaMzSg\\input\\05nc\\"
    colors=[
        "#FF0000", #Red
        "#FF7F00", #Organge
        "#FFFF00", #Yellow
        "#7FFF00", #Chartreuse
        "#00FF00", #Green
        "#00FF7F", #Spring Green
        "#00FFFF", #Cyan
        "#007FFF", #Azure
        "#0000FF", #Blue
        "#7F00FF", #Violet
        "#FF00FF", #Magenta
        "#FF007F", #Rose
    ]
    #Create a ListedColormap
    cmap=ListedColormap(colors)
    fig, axes = plt.subplots(2,2,figsize=(18,16))
    for r in range(4):
        tif=tif_dir+"MZ_season"+str(r%2+1)+"_"+tfile[r//2]+".tif"  #tif file for plotting
        with rasterio.open(tif) as src:
            raster_data = src.read(1)  # 读取第一个波段数据
            #get the bounds of the raster
            bounds=src.bounds
            #get the NoData value if it exists
            nodata_value=src.nodata
        #Mask NoData values if necessary
        if nodata_value is not None:
            raster_data=np.ma.masked_equal(raster_data,nodata_value)
        if r//2==0:
            img=axes[r//2,r%2].imshow(raster_data,cmap=cmap,vmin=1,vmax=12, extent=[bounds.left,bounds.right,bounds.bottom,
                                                                                    bounds.top])
            africa.boundary.plot(ax=axes[r//2,r%2],linewidth=0.5,color='black')
            cbar = fig.colorbar(img, ax=axes[0,r%2], ticks=np.arange(1,13),orientation='vertical', fraction=0.1, pad=0.05)
            cbar.set_label('Sowing month',fontsize=14,rotation=270,labelpad=20)
            cbar.ax.yaxis.set_label_position('right')  #Position the label on the left
            cbar.ax.yaxis.set_tick_params(labelsize=12) #Set the label size
        if r//2==1:
            #Remove pixels with value 0
            #Create a mask for values that are not zero
            #mask=raster_data!=0
            mask=np.ma.masked_where(raster_data==0,raster_data)  #create a mask where raster_data==0
            #Normalize the data
            vmin, vmax = np.percentile(mask.compressed(), (2, 98)) #clip the data  between 2nd and 98th percentiles
            norm = Normalize(vmin=vmin, vmax=vmax)
            img1=axes[r//2,r%2].imshow(raster_data,cmap="plasma",norm=norm,
                                       extent=[bounds.left,bounds.right,bounds.bottom,bounds.top])
            africa.boundary.plot(ax=axes[r//2,r%2],linewidth=0.5,color='black')
            #viridis
            cbar = fig.colorbar(img1, ax=axes[1,r%2],orientation='vertical', fraction=0.1, pad=0.05)
            cbar.set_label('Harvest area (ha)',fontsize=14,rotation=270,labelpad=20)
            cbar.ax.yaxis.set_label_position('right')  #Position the label on the left
            cbar.ax.yaxis.set_tick_params(labelsize=12) #Set the label size
        axes[r//2,r%2].set_xlabel("Longitude", fontsize=18)
        axes[r//2,r%2].set_ylabel("Latitude", fontsize=18)
        axes[r//2,r%2].set_xticks([-20,-10,0,10,20,30,40,50],
                           ['$20^\circ$W','$10^\circ$W','$0^\circ$','$10^\circ$E','$20^\circ$E','$30^\circ$',
                            '$40^\circ$E','$50^\circ$E'],fontsize=14)
        axes[r//2,r%2].set_yticks([-30,-20,-10,0,10,20,30,40],
                           ['$30^\circ$S','$20^\circ$S','$10^\circ$S','$0^\circ$','$10^\circ$N','$20^\circ$N',
                            '$30^\circ$N','$40^\circ$N'],fontsize=14)
        axes[r//2,r%2].text(-35,45,chr(97+r),fontsize=26)
    fig.savefig(fig_dir + "Supp_S1.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S2#######################################
def Supplementary_Plot_S2():
    # Define color
    # type_colors = [{-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"},{-1: "red", 0: "lightgrey", 1: "green"}] #blue, red, green, maturity type
    # change_colors = {-1: "red", 0: "lightgrey", 1: "blue"} #change
    # African country bounday
    africa = gpd.read_file("E:/RSG Dropbox/Wei Xiong/Works/CurrentProcessing/0_AfricanMaizeSorghum/Plots/CountryBoundaryLines/ne_110m_admin_0_countries.shp")
    africa = africa[africa['CONTINENT'] == 'Africa']
    tif_dir = "D:\\works\\AfricaMzSg\\input\\05nc\\"
    colors = [
        "#FF0000",  # Red
        "#FFFF00",  # Yellow
        "#FF00FF",  # Magenta  "#7FFF00", #Chartreuse
        "#00FF00",  # Green
        "#00FFFF",  # Cyan
        "#0000FF",  # Blue
    ]
    class_name = {
        1: "Wet upper mid-altitude",
        2: "Wet lower mid-altitude",
        3: "Dry mid-altitude",
        4: "Wet lowland",
        5: "Dry lowland",
        6: "Highland"
    }
    # Create a ListedColormap
    cmap = ListedColormap(colors)
    fig, ax = plt.subplots(figsize=(9, 8))
    # for r in range(4):
    tif = tif_dir + "MaizeME_Africa11.tif"  # tif file for plotting
    with rasterio.open(tif) as src:
        raster_data = src.read(1)  # 读取第一个波段数据
        # get the bounds of the raster
        # bounds=src.bounds
        # lon_range=[africa.total_bounds[0],africa.total_bounds[2]]
        # lat_range=[africa.total_bounds[1],africa.total_bounds[3]]
        bounds = src.bounds
        # show(src,ax=ax)
        # get the NoData value if it exists
        nodata_value = src.nodata
    # Mask NoData values if necessary
    if nodata_value is not None:
        raster_data = np.ma.masked_equal(raster_data, nodata_value)
    img = ax.imshow(raster_data, cmap=cmap, extent=[bounds.left, bounds.right, bounds.bottom, bounds.top])
    africa.boundary.plot(ax=ax, linewidth=0.5, color='black')
    # cbar = fig.colorbar(img, ax=ax, ticks=np.arange(1,13),orientation='vertical', fraction=0.1, pad=0.05)
    patches = [
        plt.plot([], [], marker="s", ms=10, ls="", color=colors[i], label=class_name[i + 1])[0]
        for i in range(len(colors))
    ]
    ax.legend(
        handles=patches,
        title="Maize Mega-Environments",
        bbox_to_anchor=(1.01, 1.01),
        loc='upper left',
        frameon=True
    )
    cbar.set_label('Mega-Environment', fontsize=14, rotation=270, labelpad=20)
    cbar.ax.yaxis.set_label_position('right')  # Position the label on the left
    cbar.ax.yaxis.set_tick_params(labelsize=12)  # Set the label size
    # "viridis"
    # ax.set_xlim([-17.62504269,51.13387])
    # ax.set_ylim(lat_range)
    ax.set_xlabel("Longitude", fontsize=18)
    ax.set_ylabel("Latitude", fontsize=18)
    ax.set_xticks([-20, -10, 0, 10, 20, 30, 40, 50],
                  ['$20^\circ$W', '$10^\circ$W', '$0^\circ$', '$10^\circ$E', '$20^\circ$E', '$30^\circ$',
                   '$40^\circ$E', '$50^\circ$E'], fontsize=14)
    ax.set_yticks([-30, -20, -10, 0, 10, 20, 30, 40],
                  ['$30^\circ$S', '$20^\circ$S', '$10^\circ$S', '$0^\circ$', '$10^\circ$N', '$20^\circ$N',
                   '$30^\circ$N', '$40^\circ$N'], fontsize=14)

    #plt.show()
    fig.savefig(fig_dir + "Supp_S2.png", format="png", dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S5#######################################
def create_scatter_density_plot(actual, simulated, model_name, ax, color):
    """
    创建散点密度图
    """
    # 确保数据是numpy数组
    actual = np.array(actual)
    simulated = np.array(simulated)
    # 移除NaN值
    mask = ~(np.isnan(actual) | np.isnan(simulated))
    actual = actual[mask]
    simulated = simulated[mask]
    # 计算统计指标
    rmse = np.sqrt(np.mean((actual - simulated) ** 2))
    r2 = np.corrcoef(actual, simulated)[0, 1] ** 2
    mae = np.mean(np.abs(actual - simulated))
    bias = np.mean(simulated - actual)
    # 创建散点图
    scatter = ax.scatter(actual, simulated, c=color, alpha=0.8, s=5, label=model_name)
    # 添加密度估计
    xy = np.vstack([actual, simulated])
    z = gaussian_kde(xy)(xy)
    # 按密度排序点，使高密度点显示在上层
    idx = z.argsort()
    actual, simulated, z = actual[idx], simulated[idx], z[idx]
    scatter = ax.scatter(actual, simulated, c=z, s=10, cmap='viridis', alpha=0.7)
    # 添加1:1线
    min_val = 0 #min(actual.min(), simulated.min())
    max_val = 14 #max(actual.max(), simulated.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, linewidth=2, label='1:1 line')
    # 设置坐标轴范围
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    # 添加统计信息文本框
    stats_text = f'{model_name}\nR² = {r2:.3f}\nRMSE = {rmse:.3f}\nMAE = {mae:.3f}\nBias = {bias:.3f}'
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=9)
    # 设置标签
    ax.set_xlabel('Observed Yield (t/ha)', fontsize=12)
    ax.set_ylabel('Simulated Yield (t/ha)', fontsize=12)
    ax.set_title(f'{model_name} performance', fontsize=14)
    # 添加图例
    ax.legend(loc='lower right')
    return r2, rmse, mae, bias
def create_scatter_plot(actual, simulated, model_name, ax, color):
    """
    创建散点密度图
    """
    # 确保数据是numpy数组
    actual = np.array(actual)
    simulated = np.array(simulated)
    # 移除NaN值
    mask = ~(np.isnan(actual) | np.isnan(simulated))
    actual = actual[mask]
    simulated = simulated[mask]
    # 计算统计指标
    rmse = np.sqrt(np.mean((actual - simulated) ** 2))
    r2 = np.corrcoef(actual, simulated)[0, 1] ** 2
    mae = np.mean(np.abs(actual - simulated))
    bias = np.mean(simulated - actual)
    # 创建散点图
    scatter = ax.scatter(actual, simulated, c=color, alpha=0.8, s=40, label=model_name)
    # 添加1:1线
    min_val = min(actual.min(), simulated.min())
    max_val = max(actual.max(), simulated.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, linewidth=2, label='1:1 line')
    # 设置坐标轴范围
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    # 添加统计信息文本框
    stats_text = f'{model_name}\nR² = {r2:.3f}\nRMSE = {rmse:.3f}\nMAE = {mae:.3f}\nBias = {bias:.3f}'
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=9)
    # 设置标签
    ax.set_xlabel('FAO reported Production (MT)', fontsize=12)
    ax.set_ylabel('Simulated Production (MT)', fontsize=12)
    ax.set_title(f'{model_name} performance', fontsize=14)
    # 添加图例
    ax.legend(loc='lower right')
    return r2, rmse, mae, bias
def Supplementary_Plot_S5():
    #EPIC-baseline
    folder="E:\\RSG Dropbox\\Wei Xiong\\Works\\CurrentProcessing\\0_AfricanMaizeSorghum\\"
    file=["EPIC_20crv3_obs_WithBaselineCultivar_EPIC","SimOutput_MZCER_baseline_lowerG2G3_1","SimOutput_MZIXM_baseline_lowerG2G3_1",
          "EPIC_20crv3_obs_cali_1","output_calipara_MZCER","output_calipara_MZIXM"]
    model=["EPIC","CERES-Maize","CSM-IXIM"]
    fig, axes = plt.subplots(3, 3, figsize=(18, 18))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    #comparison in grid scale
    for f in range(6):
        filename=folder+"Results\\Calibration_Step2\\"+file[f]+".csv"
        df=pd.read_csv(filename)
        if f<4:
            if filename.find("EPIC") != -1:
                df["simyld"]=df.iloc[:,9:].mean(axis=1)
                df=df[["obsy","simyld"]].dropna()
                df.columns=['obsy','simy']
            else:
                df=df[["mz_se1_yld_2007_2016","sim_yld"]].dropna()
                df.columns=['obsy','simy']
                df['simy']=df['simy']/1000
        else:
            df=df[["obsy","simy"]].dropna()
            df['simy']=df['simy']/1000
        r2_1, rmse_1, mae_1, bias_1 = create_scatter_density_plot(df['obsy'], df['simy'], model[f%3], axes[f//3,f%3], colors[0])
    #comparison in country scale
    df_fao=pd.read_csv(folder+"input\\FAOSTAT_African_Countries_data_en_1-21-2025.csv")
    yld=df_fao[(df_fao['Year']>=2007)&(df_fao['Year']<=2016)&(df_fao['Element']=='Yield')][['Area','Year','Value']]
    yld.columns=['country','year','yield']
    area=df_fao[(df_fao['Year']>=2007)&(df_fao['Year']<=2016)&(df_fao['Element']=='Area harvested')][['Area','Year','Value']]
    area.columns=['country','year','area']
    #df_fao=yld.merge(area)
    temp=area.groupby('country')['area'].max().reset_index()  #mean area
    df_fao=yld[['country','year','yield']].merge(temp)
    df_fao['prod']=(df_fao['yield']/1000)*df_fao['area']
    df_fao=df_fao.groupby('country')['prod'].mean().reset_index()
    df_fao.columns=['country','faoprod']
    #Simulated producation after yield calibration
    area=pd.read_csv(folder+"input\\Africa_SimGrid_Confirmed_5min_4calibration.csv")[['lat','lon','country','A']]
    #EPIC simulated country yield
    for f in range(6,9):
        yld=pd.read_csv(folder+"Results\\Calibration_Step2\\"+file[f-3]+".csv")
        if f==6:
            yld["simy"]=yld.iloc[:,9:].mean(axis=1)
        else:
            yld["simy"]=yld['simy']/1000
        yld=yld[['lon','lat','simy']]
        df=area.merge(yld,how='left')
        df['simprod']=df['simy']*df['A']
        df=df.groupby('country')['simprod'].sum().reset_index()
        df=df.merge(df_fao,how='left')
        df=df[df['simprod']>0].dropna()
        r2_1, rmse_1, mae_1, bias_1 = create_scatter_plot(df['faoprod']/10000000, df['simprod']/10000000, model[f%3], axes[f//3,f%3], colors[0])
    axes[0,0].text(-7.2,8,"With parameter of\n benchmark cultivars",fontsize=12)
    axes[1,0].text(-7.5,8,"With grid-calibrated\n cultivar parameters",fontsize=12)
    axes[2,0].text(-1.4,1.6,"Camparison with\n FAO production",fontsize=12)
    fig.savefig(fig_dir + "Supp_S5.png", format="png", dpi=300, bbox_inches='tight')
##########################################################FOR SUPPLEMENTARY FIGURE S6#######################################
def Supplementary_Plot_S6():
    in_dir="D:\\works\\AfricaMzSg\\input\\"
    africa=gpd.read_file("E:/RSG Dropbox/Wei Xiong/Works/CurrentProcessing/0_AfricanMaizeSorghum/Plots/CountryBoundaryLines/ne_110m_admin_0_countries.shp")
    africa=africa[africa['CONTINENT']=='Africa']
    #Plot for comparison between reported and simulated maize yields
    df0=pd.read_csv(in_dir+"mz_calibrated.csv")
    resolution = 0.5 #define the resolution of the grid
    fig, axes = plt.subplots(2,2,figsize=(18,16))
    for r in range(2):
        df=df0[(df0.obsy<15)&(df0.sea==(r+1))]
        #create a grid for the raster
        #Create the grid boundaries
        lat_min, lat_max = df['lat'].min(), df['lat'].max()
        lon_min, lon_max = df['lon'].min(), df['lon'].max()
        #Create grid points
        lat_bins = np.arange(lat_min, lat_max + resolution, resolution)
        lon_bins = np.arange(lon_min, lon_max + resolution, resolution)
        #Create a 2D histogram (rasterize the data)
        raster1, lat_edges, lon_edges = np.histogram2d(
            df['lat'], df['lon'], bins=[lat_bins, lon_bins], weights=df['obsy'])
        raster2, lat_edges, lon_edges = np.histogram2d(
            df['lat'], df['lon'], bins=[lat_bins, lon_bins], weights=df['mean'])
        #Normalize the raster by the number of points in each grid cell
        counts, _, _ = np.histogram2d(df['lat'], df['lon'], bins=[lat_bins, lon_bins])
        raster1 = np.divide(raster1, counts, out=np.full_like(raster1,np.nan), where=counts != 0)
        raster2 = np.divide(raster2, counts, out=np.full_like(raster2,np.nan), where=counts != 0)
        sc1 = axes[r,0].imshow(raster1, extent=[lon_min, lon_max, lat_min, lat_max], origin='lower', cmap='viridis', alpha=0.8)
        africa.boundary.plot(ax=axes[r,0],linewidth=0.5,color='black')
        #fig.colorbar(sc1, ax=axes[r,0], label='Observed Yields')
        axes[r,0].set_xlabel("Longitude", fontsize=18)
        axes[r,0].set_ylabel("Latitude", fontsize=18)
        axes[r,0].set_xticks([-20,-10,0,10,20,30,40,50],
                           ['$20^\circ$W','$10^\circ$W','$0^\circ$','$10^\circ$E','$20^\circ$E','$30^\circ$',
                            '$40^\circ$E','$50^\circ$E'],fontsize=14)
        axes[r,0].set_yticks([-30,-20,-10,0,10,20,30,40],
                           ['$30^\circ$S','$20^\circ$S','$10^\circ$S','$0^\circ$','$10^\circ$N','$20^\circ$N',
                            '$30^\circ$N','$40^\circ$N'],fontsize=14)
        #sc1 = axes[0].scatter(df['lon'],df['lat'], c=df['obsy'], cmap='viridis')
        #africa.boundary.plot(ax=axes[0],linewidth=0.5,color='black')
        sc2 = axes[r,1].imshow(raster2, extent=[lon_min, lon_max, lat_min, lat_max], origin='lower', cmap='viridis', alpha=0.8)
        africa.boundary.plot(ax=axes[r,1],linewidth=0.5,color='black')
        #fig.colorbar(sc2,ax=axes[r,1],label='Simulated Yields')
        axes[r,1].set_xlabel("Longitude", fontsize=18)
        axes[r,1].set_ylabel("Latitude", fontsize=18)
        axes[r,1].set_xticks([-20,-10,0,10,20,30,40,50],
                           ['$20^\circ$W','$10^\circ$W','$0^\circ$','$10^\circ$E','$20^\circ$E','$30^\circ$',
                            '$40^\circ$E','$50^\circ$E'],fontsize=14)
        axes[r,1].set_yticks([-30,-20,-10,0,10,20,30,40],
                           ['$30^\circ$S','$20^\circ$S','$10^\circ$S','$0^\circ$','$10^\circ$N','$20^\circ$N',
                            '$30^\circ$N','$40^\circ$N'],fontsize=14)
        vmin = min(np.nanmin(raster1), np.nanmin(raster2))
        vmax = max(np.nanmin(raster1), np.nanmax(raster2))
        cbar = fig.colorbar(sc1, ax=axes[r,1], orientation='vertical', fraction=0.1, pad=0.05)
        cbar.set_label('Maize Yield (t/ha)',fontsize=14,rotation=270,labelpad=20)
        cbar.ax.yaxis.set_label_position('right')  #Position the label on the left
        cbar.ax.yaxis.set_tick_params(labelsize=12) #Set the label size
    #sc2 = axes[1].scatter(df['lon'],df['lat'],c=df['mean'],cmap='viridis')
    for i in range(4):axes[i//2,i%2].text(-35,45,chr(97+i),fontsize=26)
    plt.subplots_adjust(wspace=0.03,hspace=0.16)
    #plt.tight_layout()
    fig.savefig(fig_dir+"Supp_S6.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S7#######################################
def yld2prodInAtable(fer):
    """
    A table containing area-weighted yield and production
    :param fer:
    :return: result
    """
    in_dir = "D:\\works\\AfricaMzSg\\input\\"
    # sea 1=primary 2=minor 3=total/combined
    result = pd.DataFrame(columns=['sowtype', 'type', 'cul', 'sea'] + [x for x in range(1971, 2022)])
    ###FAO Reported yield and production, computed from detrended yield, production = detrend yield  x max(area in 1971-2021)
    df = pd.read_csv(in_dir + "FAOSTAT_African_Countries_data_en_1-21-2025.csv")
    #Real FAO production of past 10 years, 2014-2023##################
    temp = df.loc[(df['Area'] == 'Africa') & (df['Year'] > 2018) & (df['Element'] == 'Production'), 'Value'].tolist()
    result.loc[len(result),] = ['fao', 'prod', 'real', '2019-2023'] + [np.nan] * (51 - len(temp)) + [x / 1000000 for x
                                                                                                     in temp]
    #########################################Estimated production by detrended FAO yield######################
    df = df[(df['Area'] == 'Africa') & (df['Year'] > 1970) & (df['Year'] < 2022)][['Element', 'Year', 'Value']]
    temp = df.loc[df['Element'] == 'Yield', ['Year', 'Value']]
    temp.columns = ['year', 'value']
    slope, intercept, r_value, p_value, std_err = stats.linregress(temp['year'], temp['value'])
    yld = ((temp['value'] + (2023 - temp['year']) * slope) / 1000).tolist()  # detrend yield to 2020 level
    result.loc[len(result),] = ['fao', 'yld', 'yld_detrend2020', ''] + [round(x, 3) for x in yld]
    prod = (np.array(yld) * df[(df['Year'] > 2011) & (df['Year'] < 2022) & (df['Element'] == 'Area harvested')][
        'Value'].mean() / 1000000).tolist()
    result.loc[len(result),] = ['fao', 'prod', 'yid_detrend2020', ''] + prod
    ####################################SIMULATED YIELD#########################################################
    # Resport sowing month and area for the primary and minor maize
    # reported sowing month
    sow = pd.read_csv(in_dir + "Africa_SimGrid_Confirmed_5min_4calibration.csv")[
        ['SIMUNIT', 'ReportedSow_se1m', 'ReportedSow_se2m', 'A']]
    sow.columns = ['SIMUNIT', 'se1m', 'se2m', 'A']
    area = pd.read_csv(in_dir + "Africa_SIMUNIT_MZ_PhysicalArea.csv")
    area = area.merge(sow[['SIMUNIT', 'A']], how='outer').fillna(0)
    area['se1a'] = area['PhysicalArea']  # Area of the primary maize is set to the physical area of maize in around 2020 (SPAM)
    area['se2a'] = area['A'] - area['PhysicalArea']  # Actual area of second maize is not know, setting to Physical-Harvest>0
    area[area < 0] = 0
    sow = sow[['SIMUNIT', 'se1m', 'se2m']]
    area = area[['SIMUNIT', 'se1a', 'se2a']]
    #######################################SIMULATED VALUE# WITHOUT FERTILIZER#################################
    sowing_window = ['ReportedSow', 'FixedSow', 'OptimumFixedSow', 'OptimumYearSow']
    in_dir = "D:/works/AfricaMzSg/output/"
    for sw in sowing_window:
        df = pd.read_csv(in_dir + "mz_yield_" + sw + "_" + fer + ".txt", delim_whitespace=True)  # change here for fertilizer
        for cul in [1, 2, 3]:
            # seson 1
            temp = df[(df['matu'] == cul) & (df['sea'] == 1)].fillna(0)
            # for season 1, we always extract yield, for season 2, we extract yield for grids with >80% yearly yield greater 1
            # temp['G1count']=[(temp.loc[i,["yld_"+str(yr) for yr in range(1971,2022)]]>=1).sum() for i in temp.index]
            # temp=temp[temp.G1count>=(51*0.8)].drop('G1count',axis=1) #only keep rows that >80% years yield greater 1
            temp = temp.merge(area, how='left').fillna(0)  # merge maize area
            yld1 = [np.dot(temp["yld_" + str(x)].values, temp['se1a'].values.T) / sum(temp['se1a']) for x in range(1971, 2022)]
            prod1 = [np.dot(temp["yld_" + str(x)].values, temp['se1a'].values.T) / 1000000 for x in range(1971, 2022)]
            result.loc[len(result),] = [sw, 'yld', cul, 1] + [round(x, 3) for x in yld1]
            result.loc[len(result),] = [sw, 'prod', cul, 1] + prod1
            a1 = sum(temp['se1a'])
            # season2
            temp = df[(df['matu'] == cul) & (df['sea'] == 2)].fillna(0)
            temp = temp.merge(area, how='left').fillna(0)
            yld2 = [np.dot(temp["yld_" + str(x)].values, temp['se2a'].values.T) / sum(temp['se2a']) for x in range(1971, 2022)]
            prod2 = [np.dot(temp["yld_" + str(x)].values, temp['se2a'].values.T) / 1000000 for x in range(1971, 2022)]
            a2 = sum(temp['se2a'])
            result.loc[len(result),] = [sw, 'yld', cul, 2] + [round(x, 3) for x in yld2]
            result.loc[len(result),] = [sw, 'prod', cul, 2] + prod2
            # aggregated yield
            yld = ((np.array(prod1) + np.array(prod2)) * 1000000 / (a1 + a2)).tolist()
            result.loc[len(result),] = [sw, 'yld', cul, 3] + [round(x, 3) for x in yld]
    return(result)
def Supplementary_Plot_S7(fer):
    type_colors = {-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"}  # blue, red, green
    result = yld2prodInAtable(fer)
    result = result[(result.sowtype != 'fao') & (result.type == 'yld') & (result.sea != 3)]
    result['mean'] = result.iloc[:, 4:].mean(axis=1)
    type = ['Planting', 'Season', 'Maturity', 'Year']
    data = []  # Create a empty list
    data.append(result[['sowtype', 'mean']].groupby('sowtype').mean().reset_index()['mean'].tolist())  # Planting
    data.append(result[['sea', 'mean']].groupby('sea').mean().reset_index()['mean'].tolist())  # Season
    data.append(result[['cul', 'mean']].groupby('cul').mean().reset_index()['mean'].tolist())  # Maturity
    data.append(result[[yr for yr in range(1971, 2022)]].mean().tolist())  # Year
    fig = plt.figure(figsize=(6, 4.5))
    plt.boxplot(data, vert=True, patch_artist=True, showmeans=False, showfliers=False,
                boxprops=dict(facecolor='lightblue', color="blue"),
                medianprops=dict(color="red", linewidth=2),
                whiskerprops=dict(color="blue", linewidth=1.5),
                capprops=dict(color="blue", linewidth=1.5))
    plt.xticks([1, 2, 3, 4], type, fontsize=12)
    plt.ylabel("Mean yield (t/ha)", fontsize=15)
    #plt.show()
    fig.savefig(fig_dir+"Supp_S7.png",format="png",dpi=300, bbox_inches='tight')
##########################################################FOR SUPPLEMENTARY FIGURE S8#######################################
def AggregatedYld(fer):
    in_dir="D:/works/AfricaMzSg/output/"
    sw=["ReportedSow","OptimumFixedSow"]
    #Maize physical area
    area=pd.read_csv("D:/works/AfricaMzSg/Input/Africa_SimGrid_Confirmed_5min_4calibration.csv")[['SIMUNIT','ReportedSow_se1m','ReportedSow_se2m','A']]
    area.columns=['SIMUNIT','se1m','se2m','A']
    yld=pd.DataFrame(columns=['sowtype','matu','sea','GT1']+[str(yr) for yr in range(1971,2022)])
    for sow in sw:
        df0=pd.read_csv(in_dir+"mz_yield_"+sow+"_"+fer+".txt",sep=r'\s+')  #Change here when you want to use wofer yield
        for cul in [1,2,3]: #three varieties
            for sea in [1,2]: #two seasons only the first season because the sowing month identifed by peaks is not feasible
                df=df0[(df0['matu']==cul)&(df0['sea']==sea)].fillna(0).reset_index()
                area1=area[['SIMUNIT','se'+str(1)+'m','A']].dropna()
                df=df.merge(area1,how='left').fillna(0)
                #not GT1
                yld.loc[len(yld)+1,]=[sow,cul,sea,0]+list((np.dot(df[['A']].values.T,
                                                    df[["yld_"+str(yr) for yr in range(1971,2022)]].values)/df['A'].sum()).flatten())
                #yes GT1 Only consider the cells whose chance of obtaining subsistent yield (>1t/ha) is greater than 80%.
                df['G1count']=[(df.loc[i,["yld_"+str(yr) for yr in range(1971,2022)]]>=1).sum() for i in df.index]
                df=df[df.G1count>=(51*0.8)] #only keep rows that >80% years yield greater 1
                yld.loc[len(yld)+1,]=[sow,cul,sea,1]+list((np.dot(df[['A']].values.T,
                                                    df[["yld_"+str(yr) for yr in range(1971,2022)]].values)/df['A'].sum()).flatten())
    return(yld)
def Supplementary_Plot_S8(fer,season):
    yld=AggregatedYld(fer)
    sw=["ReportedSow","OptimumFixedSow"]
    yld['SI']=yld[[str(yr) for yr in range(1971,2022)]].mean(axis=1)**2/(100*yld[[str(yr) for yr in range(1971,2022)]].std(axis=1))
    fig = plt.figure(figsize=(12,6))  #,width_ratios=[0.5,1]
    gs=fig.add_gridspec(2,3)
    ax1=fig.add_subplot(gs[0:2,0])
    df=yld[(yld['GT1']==0)&(yld['sea']==season)]    #Do not remove grid that over 80% years yield is greater than 1
    #SI for Africa
    data={
        'SI': df['SI'].tolist(),  #Values for the six bars
        'BarWidth': [0.7,0.7,0.7,0.7,0.7,0.7],
        'Color': ['#CCFFE7','#E3C5E7','#FFEDB7','#009951','#AF58BA','#EAAF00']
    }
    type_colors=[['#CCFFE7','#E3C5E7','#FFEDB7'],['#009951','#AF58BA','#EAAF00']] #Green, purple, pink
    type_lines=[':',"-"]
    #Green:#CCFFE7,#009951 Purple: #E3C5E7 #AF58BA Pink #FFEDB7 #EAAF00
    label=['Early','Medium','Late']
    y_p=[5,4,3,2,1,0]
    for i,(value,width,color) in enumerate(zip(data['SI'],data['BarWidth'],data['Color'])):
        ax1.barh(y_p[i],value,height=width,color=color)
        #axes[0].text(value-value/1.5,y_p[i]-0.1,label[i])
    ax1.set_yticks([0,1,2,3,4,5],['Late','Medium','Early','Late','Medium','Early'])
    ax1.set_ylabel('Optimum Sow            Reported Sow',fontsize=15)
    ax1.set_xlabel('Suitability Index (SI)',fontsize=15)
    #For decades
    ax2=fig.add_subplot(gs[0,1:3])
    sow=0
    for cul in range(3):
        #SI Change  - yld
        y = df.loc[(df.matu==cul+1)&(df.sowtype==sw[sow])&(df.sea==season),[str(x) for x in range(1971,2022)]]
        y = pd.DataFrame(y.values[0])
        if sow==1: y=y-4
        #plt.plot()
        y = (y[0].groupby(np.arange(len(y)) // 10).mean())**2/(100*y[0].groupby(np.arange(len(y)) // 10).std())
        y = y.dropna()
        x = list(range(1975,2021,10))
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        ax2.plot(x,y,marker='o',markersize=10,color=type_colors[1][cul],linestyle=":",label=label[cul]) #,
        #trend
        z = np.polyfit(x,y,1) #fit a 1st-degree polynomial (linear)
        p = np.poly1d(z)
        ax2.plot(x, p(x), "-", color=type_colors[1][cul])
        ax2.set_xticks([1976,1986,1996,2006,2016])
        ax2.set_xticklabels([])
        ax2.set_ylabel("SI",fontsize=15,labelpad=10)
    ax3=fig.add_subplot(gs[1,1:3])
    sow=1
    for cul in range(3):
        #SI Change  - yld
        y = df.loc[(df.matu==cul+1)&(df.sowtype==sw[sow])&(df.sea==season),[str(x) for x in range(1971,2022)]]
        y = pd.DataFrame(y.values[0])
        #plt.plot()
        y = (y[0].groupby(np.arange(len(y)) // 10).mean())**2/(100*y[0].groupby(np.arange(len(y)) // 10).std())
        y = y.dropna()
        x = list(range(1975,2021,10))
        ax3.plot(x,y,marker='o',markersize=10,color=type_colors[sow][cul],linestyle=":") #,
        #trend
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        z = np.polyfit(x,y,1) #fit a 1st-degree polynomial (linear)
        p = np.poly1d(z)
        ax3.plot(x, p(x), "-", color=type_colors[sow][cul])
        ax3.set_xticks([1976,1986,1996,2006,2016])
        ax3.set_xticklabels(["1970s","1980s","1990s","2000s","2010s"])
        ax3.set_xlabel("Decade",fontsize=15)
        ax3.set_ylabel("SI",fontsize=15)
        ax3.tick_params(axis='both', which='major')
    ax2.legend(bbox_to_anchor=(1.2, 1.03), loc='upper right',fontsize=10)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.30, hspace=0.15)
    ax1.text(-0.6,5.8,'a',fontsize=20)
    ax2.text(1968,ax2.get_ylim()[1]+(ax2.get_ylim()[1]-ax2.get_ylim()[0])/20,'b',fontsize=20)
    ax3.text(1968,ax3.get_ylim()[1]+(ax3.get_ylim()[1]-ax3.get_ylim()[0])/20,'c',fontsize=20)
    fig.savefig(fig_dir+"\\Supp_S8.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)

##########################################################FOR SUPPLEMENTARY FIGURE S9#######################################
def Supplementary_Plot_S9(crop,fer,G1):
    point=[[151215226.6,2.708938],[48290267.12,5.104918]]  #season 1 neib=24, nmarket=31, season 2 neib=10, nmarket=4
    in_dir="D:/works/AfricaMzSg/results/"+crop+"/"+fer+"/"
    df=pd.read_csv(in_dir+crop+"_ProdBenefitsToSegmentationPara_"+fer+G1+".csv")
    fig, ax = plt.subplots(1, 2, figsize=(12,4))
    for se in range(2):
        #base=df.loc[(df.season==(se+1))&(df.nneigh==0)&(df.nseg==0)]
        temp=df.loc[(df.season==(se+1))&(df.nneigh>0)&(df.nseg>0)]
        #temp['ch_prod']=100*(temp['prod_mean']-base.iloc[0,3])/base.iloc[0,3]
        #temp['ch_cv']=100*(temp['prod_cv']-base.iloc[0,4])/base.iloc[0,4]
        sc=ax[se].scatter(temp['prod_mean'],temp['prod_cv'], c=temp['nseg'],cmap='viridis',s=10,alpha=0.9)
        ax[se].scatter(point[se][0],point[se][1],marker='s',s=60,edgecolor='red',facecolor='none',linewidth=2)
        ax[se].set_ylabel("Production CV (%)",fontsize=12)
        ax[se].set_xlabel("Production mean (T)", fontsize=12)
        ax[se].text(ax[se].get_xlim()[0]-(ax[se].get_xlim()[1]-ax[se].get_xlim()[0])/7,
            ax[se].get_ylim()[1]+(ax[se].get_ylim()[1]-ax[se].get_ylim()[0])/20,chr(97+se),fontsize=18)
        #ax[se].axvline(x=0,color='k',linestyle='--',linewidth=0.5,alpha=0.8)
        #ax[se].axhline(y=0,color='k',linestyle='--',linewidth=0.5,alpha=0.8)
    #fig.subplots_adjust(right=0.85,wspace=0.3)
    ax[0].text(106000000,2.9,"n_neighbors = "+str(para[0][0])+"\nn_clusters = "+str(para[0][1]),fontsize=12,style='italic',c='red')
    ax[1].text(39500000,5.0,"n_neighbors = "+str(para[1][0])+"\nn_clusters = "+str(para[1][1]),fontsize=12,style='italic',c='red')
    
    cbar_ax=fig.add_axes([0.91,0.110,0.016,0.77])
    cbar=plt.colorbar(sc,cax=cbar_ax)
    cbar.set_label('Number of market',fontsize=12,rotation=270,labelpad=15)
    fig.savefig(fig_dir+"\\Supp_S9.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S10#######################################
def Supplementary_Plot_S10(fer,season):
    yld=AggregatedYld(fer)
    sw=["ReportedSow","OptimumFixedSow"]
    yld['SI']=yld[[str(yr) for yr in range(1971,2022)]].mean(axis=1)**2/(100*yld[[str(yr) for yr in range(1971,2022)]].std(axis=1))
    fig = plt.figure(figsize=(12,6))  #,width_ratios=[0.5,1]
    gs=fig.add_gridspec(2,3)
    ax1=fig.add_subplot(gs[0:2,0])
    df=yld[(yld['GT1']==0)&(yld['sea']==season)]    #Do not remove grid that over 80% years yield is greater than 1
    #SI for Africa
    data={
        'SI': df['SI'].tolist(),  #Values for the six bars 
        'BarWidth': [0.7,0.7,0.7,0.7,0.7,0.7],
        'Color': ['#CCFFE7','#E3C5E7','#FFEDB7','#009951','#AF58BA','#EAAF00']
    }
    type_colors=[['#CCFFE7','#E3C5E7','#FFEDB7'],['#009951','#AF58BA','#EAAF00']] #Green, purple, pink
    type_lines=[':',"-"]
    #Green:#CCFFE7,#009951 Purple: #E3C5E7 #AF58BA Pink #FFEDB7 #EAAF00
    label=['Early','Medium','Late']
    y_p=[5,4,3,2,1,0]
    for i,(value,width,color) in enumerate(zip(data['SI'],data['BarWidth'],data['Color'])):
        ax1.barh(y_p[i],value,height=width,color=color)
        #axes[0].text(value-value/1.5,y_p[i]-0.1,label[i])
    ax1.set_yticks([0,1,2,3,4,5],['Late','Medium','Early','Late','Medium','Early'])
    ax1.set_ylabel('Optimum Sow            Reported Sow',fontsize=15)
    ax1.set_xlabel('Suitability Index (SI)',fontsize=15)
    #For decades
    ax2=fig.add_subplot(gs[0,1:3])
    sow=0
    for cul in range(3):
        #SI Change  - yld
        y = df.loc[(df.matu==cul+1)&(df.sowtype==sw[sow])&(df.sea==season),[str(x) for x in range(1971,2022)]]
        y = pd.DataFrame(y.values[0])
        if sow==1: y=y-4
        #plt.plot()
        y = (y[0].groupby(np.arange(len(y)) // 10).mean())**2/(100*y[0].groupby(np.arange(len(y)) // 10).std())
        y = y.dropna()
        x = list(range(1975,2021,10))
        ax2.plot(x,y,marker='o',markersize=10,color=type_colors[1][cul],linestyle=":",label=label[cul]) #,
        #trend
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        z = np.polyfit(x,y,1) #fit a 1st-degree polynomial (linear)
        p = np.poly1d(z)
        ax2.plot(x, p(x), "-", color=type_colors[1][cul])
        ax2.set_xticks([1976,1986,1996,2006,2016])
        ax2.set_xticklabels([])
        ax2.set_ylabel("SI",fontsize=15,labelpad=10)
    ax3=fig.add_subplot(gs[1,1:3])
    sow=1
    for cul in range(3):
        #SI Change  - yld
        y = df.loc[(df.matu==cul+1)&(df.sowtype==sw[sow])&(df.sea==season),[str(x) for x in range(1971,2022)]]
        y = pd.DataFrame(y.values[0])
        #plt.plot()
        y = (y[0].groupby(np.arange(len(y)) // 10).mean())**2/(100*y[0].groupby(np.arange(len(y)) // 10).std())
        y = y.dropna()
        x = list(range(1975,2021,10))
        ax3.plot(x,y,marker='o',markersize=10,color=type_colors[sow][cul],linestyle=":") #,
        #trend
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        z = np.polyfit(x,y,1) #fit a 1st-degree polynomial (linear)
        p = np.poly1d(z)
        ax3.plot(x, p(x), "-", color=type_colors[sow][cul])
        ax3.set_xticks([1976,1986,1996,2006,2016])
        ax3.set_xticklabels(["1970s","1980s","1990s","2000s","2010s"])
        ax3.set_xlabel("Decade",fontsize=15)
        ax3.set_ylabel("SI",fontsize=15)
        ax3.tick_params(axis='both', which='major')
    ax2.legend(bbox_to_anchor=(1.2, 1.03), loc='upper right',fontsize=10)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.30, hspace=0.15)
    ax1.text(-0.6,5.8,'a',fontsize=20)
    ax2.text(1968,ax2.get_ylim()[1]+(ax2.get_ylim()[1]-ax2.get_ylim()[0])/20,'b',fontsize=20)
    ax3.text(1968,ax3.get_ylim()[1]+(ax3.get_ylim()[1]-ax3.get_ylim()[0])/20,'c',fontsize=20)
    fig.savefig(fig_dir+"\\Supp_Fig_S10.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S11#######################################
def Supplementary_Plot_S11(crop,fer,G1):
    type_colors = ['#009951', '#AF58BA', '#EAAF00']
    label = ['Early', 'Medium', 'Late']
    cul = ['.1', '0', '1']
    fig, ax = plt.subplots(3, 1, figsize=(8, 10))
    # Tropical Bimodel
    df0 = df.loc[df['SIMUNIT'] == 21599]
    text_pos = [[9, 1.5], [6.9, 1.0], [10.9, 1.3]]
    for c in range(3):
        # season 1
        df1 = df0[[col for col in df0.columns if 'm' + cul[c] in col]].mean(axis=0).tolist()
        ax[0].plot(range(1, 13), df1, color=type_colors[c], label=label[c])
        df1 = df1[:8]
        mxvalue = max(df1)
        mon = df1.index(mxvalue) + 1
        ax[0].plot([mon, mon], [-0.2, mxvalue], linestyle='dashed', color='black', linewidth=0.5)
        # if mon==12: mon=1
        # ax.plot([mon,mon+3+c],[0+c*0.2]*2,color=type_colors[c],linewidth=10,alpha=0.5)
        df2 = df0['m' + cul[c] + 's' + str(mon).zfill(2)].values
        si = round(df2.mean() ** 2 / (100 * df2.std()), 4)
        ax[0].text(mon, 5 - mxvalue, "SI=" + str(si), color=type_colors[c])
    c = 1
    ax[0].plot([mon + 0.30, mon + 3 + c - 0.25], [0.2 + c * 0.2] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[0].text(3.5, 0.3, "Primary season")
    # season 2
    for c in range(3):
        # list(set([i for i in range(1,13)])-set([i for i in range(mon,mon+5)])) #other months expect primary season
        df1 = df0[[col for col in df0.columns if 'm' + cul[c] in col]].mean(axis=0).tolist()
        sea1_mon = [i - 1 for i in range(mon, mon + 3 + c)]  # season 1 month
        df1 = [0 if i in sea1_mon else df1[i] for i in range(12)]
        mxvalue1 = max(df1)
        mon1 = df1.index(mxvalue1) + 1
        ax[0].plot([mon1, mon1], [-0.2, mxvalue1], linestyle='dashed', color='black', linewidth=0.5)
        # if mon==12: mon=1
        # ax.plot([mon,mon+3+c],[0+c*0.2]*2,color=type_colors[c],linewidth=10,alpha=0.5)
        df2 = df0['m' + cul[c] + 's' + str(mon1).zfill(2)].values
        si = round(df2.mean() ** 2 / (100 * df2.std()), 4)
        ax[0].text(mon1, mxvalue1 - c, "SI=" + str(si), color=type_colors[c])
    c = 2
    ax[0].plot([8 + 0.25, 8 + 3 + 1 - 0.25], [0 + 2 * 0.2] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[0].plot([1.25, 2 - 0.25], [0 + 2 * 0.2] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[0].text(9.2, 0.3, "Minor season")
    ax[0].set_ylim([-0.2, 4.3])
    ax[0].set_xticks([i for i in range(1, 13)], [i for i in range(1, 13)])
    ax[0].set_yticks([i for i in range(5)], [i for i in range(5)])
    # ax[0].set_xlabel("Month",fontsize=12)
    ax[0].set_ylabel("Yield (t/ha)", fontsize=12)
    ax[0].legend(bbox_to_anchor=(1.19, 1.02), loc='upper right', fontsize=10)
    ##Temporal Unimodel low rainfall
    df0 = df.loc[df['SIMUNIT'] == 58260]
    for c in range(3):
        # season 1
        df1 = df0[[col for col in df0.columns if 'm' + cul[c] in col]].mean(axis=0).tolist()
        ax[1].plot(range(1, 13), df1, color=type_colors[c], label=label[c])
        mxvalue = max(df1)
        mon = df1.index(mxvalue) + 1
        ax[1].plot([mon, mon], [-0.2, mxvalue], linestyle='dashed', color='black', linewidth=0.5)
        # if mon==12: mon=1
        # ax.plot([mon,mon+3+c],[0+c*0.2]*2,color=type_colors[c],linewidth=10,alpha=0.5)
        df2 = df0['m' + cul[c] + 's' + str(mon).zfill(2)].values
        si = round(df2.mean() ** 2 / (100 * df2.std()), 4)
        ax[1].text(mon, mxvalue / 0.5 - 4, "SI=" + str(si), color=type_colors[c])
    c = 0
    ax[1].plot([mon + 0.30, mon + 3 + c - 0.25], [0.2 + c * 0.2] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[1].text(6.5, 0.15, "Primary season")
    ax[1].set_ylim([-0.2, 3.2])
    ax[1].set_xticks([i for i in range(1, 13)], [i for i in range(1, 13)])
    ax[1].set_yticks([i for i in range(4)], [i for i in range(4)])
    # ax[1].set_xlabel("Month",fontsize=12)
    ax[1].set_ylabel("Yield (t/ha)", fontsize=12)
    ##Temporal Unimodel high rainfall with long period
    df0 = df.loc[df['SIMUNIT'] == 162680]
    for c in range(3):
        # season 1
        df1 = df0[[col for col in df0.columns if 'm' + cul[c] in col]].mean(axis=0).tolist()
        ax[2].plot(range(1, 13), df1, color=type_colors[c], label=label[c])
        mxvalue = max(df1)
        mon = df1.index(mxvalue) + 1
        ax[2].plot([mon, mon], [-0.2, mxvalue], linestyle='dashed', color='black', linewidth=0.5)
        # if mon==12: mon=1
        # ax.plot([mon,mon+3+c],[0+c*0.2]*2,color=type_colors[c],linewidth=10,alpha=0.5)
        df2 = df0['m' + cul[c] + 's' + str(mon).zfill(2)].values
        si = round(df2.mean() ** 2 / (100 * df2.std()), 4)
        if c == 1:
            ax[2].text(mon - 1, mxvalue - 0.3, "SI=" + str(si), color=type_colors[c])
        else:
            ax[2].text(mon, mxvalue - 0.2, "SI=" + str(si), color=type_colors[c])
    c = 2
    ax[2].plot([mon + 0.25, mon + 1 - 0.25], [0.3] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[2].plot([1 + 0.25, 4 + 1 - 0.25], [0.3] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[2].text(1.7, 0.25, "Primary season")
    # season 2
    for c in range(3):
        # list(set([i for i in range(1,13)])-set([i for i in range(mon,mon+5)])) #other months expect primary season
        df1 = df0[[col for col in df0.columns if 'm' + cul[c] in col]].mean(axis=0).tolist()
        sea1_mon = [i - 1 for i in range(mon, mon + 3 + c)]  # season 1 month
        sea1_mon = [x % 12 for x in sea1_mon] + [i for i in range(mon - 3 - c, mon + 1)]
        df1 = [0 if i in sea1_mon else df1[i] for i in range(12)]
        mxvalue1 = max(df1)
        mon1 = df1.index(mxvalue1) + 1
        ax[2].plot([mon1, mon1], [-0.2, mxvalue1], linestyle='dashed', color='black', linewidth=0.5)
        # if mon==12: mon=1
        # ax.plot([mon,mon+3+c],[0+c*0.2]*2,color=type_colors[c],linewidth=10,alpha=0.5)
        df2 = df0['m' + cul[c] + 's' + str(mon1).zfill(2)].values
        si = round(df2.mean() ** 2 / (100 * df2.std()), 4)
        ax[2].text(mon1, mxvalue1 - c, "SI=" + str(si), color=type_colors[c])
    c = 2
    ax[2].plot([mon1 + 0.25, mon1 + c + 3 - 0.25], [0.2] * 2, color=type_colors[c], linewidth=20, alpha=0.5)
    ax[2].text(7.8, 0.15, "Minor season")
    ax[2].set_ylim([-0.2, 5.2])
    ax[2].set_xticks([i for i in range(1, 13)], [i for i in range(1, 13)])
    ax[2].set_yticks([i for i in range(6)], [i for i in range(6)])
    ax[2].set_xlabel("Month", fontsize=12)
    ax[2].set_ylabel("Yield (t/ha)", fontsize=12)
    ax[0].text(-0.5, 4.5, "a", fontsize=15)
    ax[1].text(-0.5, 3.3, "b", fontsize=15)
    ax[2].text(-0.5, 5.3, "c", fontsize=15)
    fig.savefig(fig_dir + "Supp_S11.png", format="png", dpi=300, bbox_inches='tight', pad_inches=0)
##########################################################FOR SUPPLEMENTARY FIGURE S12#######################################
def Supplementary_Plot_S12(crop,fer,G1):
    cmap1 = plt.cm.tab20
    colors=cmap1.colors
    fig, ax = plt.subplots(2, 1, figsize=(12,12))
    in_dir="D:/works/AfricaMzSg/results/"+crop+"/"+fer+"/"
    fig_dir="E:/RSG Dropbox/Wei Xiong/Works/CurrentProcessing/0_AfricanMaizeSorghum/NF_NATFOOD-20501074/R1/Plots/"
    seg=pd.read_csv(in_dir+"mz_CountryProdBenefitByClustering_"+fer+G1+".csv") 
    for se in [1,3]:
        df1=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wcluster'),]
        df2=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wocluster'),]
        df2=df1[['country']].merge(df2,on='country',how='left')
        df=pd.DataFrame(df1.iloc[:,4:55].values-df2.iloc[:,4:55].values)/1000000
        df1=df1[['country','nmarket']]
        df1['benefit']=df.mean(axis=1).tolist()
        df1['benefit_std']=df.std(axis=1).tolist()
        data=df1[df1.benefit!=0].sort_values(by='benefit',ascending=False).replace("United Republic of Tanzania",
                                                   "Tanzania").replace("Democratic Republic of the Congo", 
                                                                       "DRC").replace("Libyan Arab Jamahiriya",
                                                                                      "Libya").replace("C??te d'Ivoire","Cote d'Ivoire")        
        bars=ax[se//2].bar(data['country'],data['benefit'],
                           yerr=data['benefit_std'],capsize=2,error_kw={'elinewidth':1,'alpha':0.6,},
                           edgecolor='black',color=[cmap1(n) for n in data['nmarket']],linewidth=0.5,width=1)
        ax[se//2].set_xticks(data['country'], data['country'], rotation=270, ha='right',fontsize=9)
        ax[se//2].set_ylabel('Production benefit (MT)',fontsize=15)
        nmarket=data['nmarket'].tolist()
        nmarket=list(set(nmarket))
        cbar_ax = fig.add_axes([0.91,0.565-(se//2)*0.455,0.02,0.315])
        nmarket_cmap=mcolors.ListedColormap([cmap1(n) for n in nmarket])
        boundaries=np.arange(0,len(nmarket)+1)
        norm=mcolors.BoundaryNorm(boundaries,len(nmarket),clip=True)
        sm=plt.cm.ScalarMappable(cmap=nmarket_cmap,norm=norm)
        cbar=fig.colorbar(sm,orientation="vertical",cax=cbar_ax)
        cbar.set_ticks([0.5+x for x in range(len(nmarket))])
        cbar.set_ticklabels([f'{i}' for i in nmarket])
        cbar.set_label("Number of markets",fontsize=14,rotation=270,labelpad=15)
        cbar.ax.tick_params(which='minor', color='white', labelcolor='white')
        ax[se//2].text(ax[se//2].get_xlim()[0]-(ax[se//2].get_xlim()[1]-ax[se//2].get_xlim()[0])/15,
                ax[se//2].get_ylim()[1]+(ax[se//2].get_ylim()[1]-ax[se//2].get_ylim()[0])/20,chr(97+se//2),fontsize=18)
    #production benefit for the minor season comes from the potential area expansion
    ax_inset=ax[1].inset_axes([0.55,0.35,0.35,0.55])
    data=pd.DataFrame(columns=['whicharea','benefit','benefit_std'])
    area=['Current','Potential']
    for se in [2,3]:
        df1=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wcluster'),]
        df2=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wocluster'),]
        df2=df1[['country']].merge(df2,on='country',how='left')
        df=pd.DataFrame(df1.iloc[:,4:55].values-df2.iloc[:,4:55].values)/1000000
        data.loc[len(data),]=[area[se-2],df.mean(axis=1).sum(),df.sum(axis=0).std()]
    ax_inset.bar(data.whicharea.tolist(),data.benefit.values,yerr=data.benefit_std.values,capsize=5,
                 width=0.6,color=['skyblue','lightblue'],edgecolor='black',linewidth=0.5)
    ax_inset.set_ylabel('Production benefit (MT)',fontsize=10)
    ax_inset.set_xticks([42,43],data.whicharea.tolist(),fontsize=10)
    ax_inset.set_xlabel('Possibile area for the minor season',fontsize=10)
    ax_inset.text(41.85,4.6,str(round(data.iloc[0,1],3))+'MT',fontsize=10)
    ax_inset.text(42.85,18,str(round(data.iloc[1,1],3))+'MT',fontsize=10)
    plt.subplots_adjust(hspace=0.45)
    fig.savefig(fig_dir+"\\Supp_S12.png",format="png",dpi=300,bbox_inches='tight', pad_inches=0)

##########################################################FOR SUPPLEMENTARY FIGURE S14#######################################
def yld2prodInAtable(crop,fer,cropmodel):
    in_dir = "D:\\works\\AfricaMzSg\\input\\"
    result = pd.DataFrame(columns=['sowtype', 'type', 'cul', 'sea'] + [x for x in range(1971, 2022)])
    df = pd.read_csv(in_dir + "FAOSTAT_African_Countries_data_en_1-21-2025.csv")
    temp = df.loc[(df['Area'] == 'Africa') & (df['Year'] > 2018) & (df['Element'] == 'Production'), 'Value'].tolist()
    result.loc[len(result),] = ['fao', 'prod', 'real', '2019-2023'] + [np.nan] * (51 - len(temp)) + [x / 1000000 for x in temp]
    #########################################Estimated production by detrended FAO yield######################
    df = df[(df['Area'] == 'Africa') & (df['Year'] > 1970) & (df['Year'] < 2022)][['Element', 'Year', 'Value']]
    temp = df.loc[df['Element'] == 'Yield', ['Year', 'Value']]
    temp.columns = ['year', 'value']
    slope, intercept, r_value, p_value, std_err = stats.linregress(temp['year'], temp['value'])
    yld = ((temp['value'] + (2023 - temp['year']) * slope) / 1000).tolist()  # detrend yield to 2020 level
    result.loc[len(result),] = ['fao', 'yld', 'yld_detrend2020', ''] + [round(x, 3) for x in yld]
    prod = (np.array(yld) * df[(df['Year'] > 2011) & (df['Year'] < 2022) & (df['Element'] == 'Area harvested')][
        'Value'].mean() / 1000000).tolist()
    result.loc[len(result),] = ['fao', 'prod', 'yid_detrend2020', ''] + prod
    ####################################SIMULATED YIELD#########################################################
    # Resport sowing month and area for the primary and minor maize reported sowing month
    sow = pd.read_csv(in_dir + "Africa_SimGrid_Confirmed_5min_4calibration.csv")[
        ['SIMUNIT', 'ReportedSow_se1m', 'ReportedSow_se2m', 'A']]
    sow.columns = ['SIMUNIT', 'se1m', 'se2m', 'A']
    area = pd.read_csv(in_dir + "Africa_SIMUNIT_MZ_PhysicalArea.csv")
    area = area.merge(sow[['SIMUNIT', 'A']], how='outer').fillna(0)
    area['se1a'] = area['PhysicalArea']  # Area of the primary maize is set to the physical area of maize in around 2020 (SPAM)
    area['se2a'] = area['A'] - area['PhysicalArea']  # Actual area of second maize is not know, setting to Physical-Harvest>0
    area[area < 0] = 0
    sow = sow[['SIMUNIT', 'se1m', 'se2m']]
    area = area[['SIMUNIT', 'se1a', 'se2a']]
    #######################################SIMULATED VALUE# WITHOUT FERTILIZER#################################
    sowing_window = ['ReportedSow', 'FixedSow', 'OptimumFixedSow', 'OptimumYearSow']
    in_dir = "D:/works/AfricaMzSg/output/"
    for sw in sowing_window:
        df = pd.read_csv(in_dir +cropmodel+"_"+crop+"_yield_" + sw + "_" + fer + ".txt", sep=r'\s+')  # change here for fertilizer
        for cul in [1, 2, 3]:
            # seson 1
            temp = df[(df['matu'] == cul) & (df['sea'] == 1)].fillna(0)
            # for season 1, we always extract yield, for season 2, we extract yield for grids with >80% yearly yield greater 1
            # temp['G1count']=[(temp.loc[i,["yld_"+str(yr) for yr in range(1971,2022)]]>=1).sum() for i in temp.index]
            # temp=temp[temp.G1count>=(51*0.8)].drop('G1count',axis=1) #only keep rows that >80% years yield greater 1
            temp = temp.merge(area, how='left').fillna(0)  # merge maize area
            prod1 = [np.dot(temp["yld_" + str(x)].values, temp['se1a'].values.T) / 1000000 for x in range(1971, 2022)]
            result.loc[len(result),] = [sw, 'prod', cul, 1] + prod1
            a1 = sum(temp['se1a'])
    return(result)

def Supplementary_Plot_S14(crop,fer):
    fig_dir=
    models=["EPIC","MZCER","MZIXM"]
    result=yld2prodInAtable(crop,fer,models[0])
    result['model']=models[0]
    for m in models[1:]:
        temp=yld2prodInAtable(crop,fer,m)
        temp=temp.iloc[3:,:]
        temp['model']=m
        result=pd.concat([result,temp])
    sw=['ReportedSow','FixedSow','OptimumYearSow'] 
    type_colors = ["#00CD6C", "#AF58BA","#FFC61E"] #blue, red, green
    xtick = [1,3,4,5,7,8,9,11,12,13]
    bar_width=0.66
    fig, ax = plt.subplots(figsize=(8,6)) 
    #fao
    mean=np.mean(result[(result['sowtype']=="fao")&(result['type']=="prod")&(result['cul']=="real")][[2017,2018,2019,2020,2021]].values)
    sd=np.std(result[(result['sowtype']=="fao")&(result['type']=="prod")&(result['cul']=="real")][[2017,2018,2019,2020,2021]].values)
    ax.bar(xtick[0],mean,bar_width,yerr=sd,color='grey',error_kw=dict(lw=1, capsize=2, capthick=1))
    #simulated
    for model in range(3):
        for ty in range(3): #sowing type
            for cul in range(3): #variety
                mean=np.mean(result.loc[(result['model']==models[model])&(result['sowtype']==sw[ty])&(result['cul']==(cul+1)),result.columns[4:54]].values)
                sd=np.std(result.loc[(result['model']==models[model])&(result['sowtype']==sw[ty])&(result['cul']==(cul+1)),result.columns[4:54]].values)
                #axes[1].bar(xtick[ty]-bar_width/3+(bar_width/3)*(cul-1),mean2.mean(),bar_width/3,color=type_colors[cul-1])
                if model<1:
                    ax.bar(xtick[1+model*3+ty]-bar_width/3+(bar_width/3)*cul,mean,bar_width/3,color=type_colors[cul],yerr=sd,
                           error_kw=dict(lw=1, capsize=2, capthick=1))
                elif model==1:
                    ax.bar(xtick[1+model*3+ty]-bar_width/3+(bar_width/3)*cul,(mean/1000)*0.6,bar_width/3,color=type_colors[cul],yerr=(sd/1000)*0.6,
                           error_kw=dict(lw=1, capsize=2, capthick=1))
                else:
                    ax.bar(xtick[1+model*3+ty]-bar_width/3+(bar_width/3)*cul,mean/1000,bar_width/3,color=type_colors[cul],yerr=sd/1000,
                           error_kw=dict(lw=1, capsize=2, capthick=1))
    for x in [2,6,10]: ax.axvline(x=x,color='k',linestyle='--',linewidth=0.5,alpha=0.8)
    ax.set_xticks([1,3,4,5,7,8,9,11,12,13],['Reported FAO',
                                            'Reported sow','Optimum sow','Varied optimum sow',
                                            'Reported sow','Optimum sow','Varied optimum sow',
                                            'Reported sow','Optimum sow','Varied optimum sow'],rotation=270,fontsize=9)
    ax.text(0.6,350,'FAO',fontsize=10)
    ax.text(3.5,350,'EPIC',fontsize=10)
    ax.text(7,350,'CERES-MAIZE',fontsize=10)
    ax.text(11,350,'CSM-IXIM',fontsize=10)
    ax.set_ylabel('Production benefit (MT)',fontsize=15)
    legend=[Patch(facecolor=type_colors[0],edgecolor='k',linewidth=1,label='Early'),
           Patch(facecolor=type_colors[1],edgecolor='k',linewidth=1,label='Medium'),
           Patch(facecolor=type_colors[2],edgecolor='k',linewidth=1,label='Late')]
    ax.legend(handles=legend,bbox_to_anchor=(1.2, 1.015), loc='upper right',fontsize=10)
    fig.savefig(fig_dir+"\\Supp_S14.png",format="png",dpi=300,bbox_inches='tight', pad_inches=0)

#Main function
fig_dir="E:/RSG Dropbox/Wei Xiong/Works/CurrentProcessing/0_AfricanMaizeSorghum/NF_NATFOOD-20501074/R1/Plots//"
crop='mz'
fer="wfer_gridcalibratedWaHi"
G1="_G1"
#Plot1 in Supplementary:
Supplementary_Plot_S1()  #Reported sowing month and harvest area
Supplementary_Plot_S2()  #Maize Mega-environment in Africa
Supplementary_Plot_S5()  #Grid simulation performance with benchmark varieties vs. grid-calibrated para.
Supplementary_Plot_S6()  #Spatial comparison between reported and simulated yields
Supplementary_Plot_S7(fer) #Yield variation caused by different factors
Supplementary_Plot_S8(fer，season=2) #SI and decadal changes for the minor season
Supplementary_Plot_S9(crop,fer,0.2) #sensitivity of production to clustering parameters
Supplementary_Plot_S10("wofer_gridcalibratedWaHi",1) #SI and decadal changes for primary season without low fertilizer input
Supplementary_Plot_S11(crop,"wofer_gridcalibratedWaHi") #cultivar distribution with low-input fertilizer scenarios
Supplementary_Plot_S12(fer) #Approach to identify primary and minor season
Supplementary_Plot_S12(crop,fer,G1)  #benefits across countries
Supplementary_Plot_S14(fer) #Production benefits estimated across crop models