"""
step 1: segmentation seed by maturity group
step 2: combine the optimum sowing month
plot  -  a. variety maturity market segmentation
         b. the optimum sowing month
         c. how much segmentation for the top countries 
         d. production benifits for the top counties
        
Plot 4 - Cultivar segmentation and its benefits to maize production
         4a. left: potential benefits by optimizing seed
             right:relationsip between number of segmented market and changes in prodution (Season 1)
         4b. the best segmentation scenario
         4c. which country has the highest benefits from market segmentation
"""
import matplotlib as mpl
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import scipy.stats as stats
import os
import seaborn as sns
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
import itertools
from statistics import mode
import squarify
import matplotlib.colors as mcolors
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"E:\RSG Dropbox\Wei Xiong\Works\CurrentProcessing\0_AfricaMaizeSorghum\Data")  #SCRIPT_DIR.parent / "data"
RESULTS_ROOT = Path(r"D:\Works\AfricaMzSg\results")   #SCRIPT_DIR.parent / "results"
SIMOUT_DIR = Path(r"D:\Works\AfricaMzSg\simout") #SCRIPT_DIR.parent / "simout"

BOUNDARY_DIR = Path(r"E:\RSG Dropbox\Wei Xiong\Works\CurrentProcessing\0_AfricanMaizeSorghum\Plots")
BOUNDARY_FILE = BOUNDARY_DIR / "CountryBoundaryLines" / "ne_110m_admin_0_countries.shp"

"""
Step 1 - Estimate production benefits for different number of market
"""
def clustering(df,nneigh,nmarket):
    # 假设列: lon, lat, variety, month
    # 2. 特征预处理
    # 提取属性特征进行聚类
    X_attr = df[['c', 'm']].values
    # 提取空间坐标
    coords = df[['x', 'y']].values
    # 标准化属性特征（防止月份的1-12比品种的-1-1权重过大）
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_attr)
    # 3. 创建空间连通性约束 (这是保证“空间相连”的关键)
    # 我们定义每个点只与其最近的 6 个邻居点连通（模拟蜂窝状或网格状布局）
    connectivity = kneighbors_graph(coords, n_neighbors=nneigh, include_self=False)
    # 4. 执行带约束的层级聚类
    # n_clusters: 你想要分割出的市场个数
    # linkage: 'ward' 最小化方差，使每个区域内部属性最接近
    model = AgglomerativeClustering(
        n_clusters=nmarket, 
        connectivity=connectivity, 
        linkage='ward'
    )
    df['seg'] = model.fit_predict(X_scaled)
    #using the mode of c and m for each segmentation
    temp=df[['seg','c','m']].groupby(['seg']).agg(mode).reset_index()  #segmentation 
    df=df[['x','y','seg']].merge(temp,on='seg',how='left')  #merge c and m for seg
    return(df)

def IdentifyCulGroup(group):
    c=group['c'].tolist()[0]  #c=str(int(group['c'][0])).replace("-",".")
    m=group['m'].tolist()[0]  #m=str(int(group['m'][0])).zfill(2)
    a=group['A'].mean()
    if m=='' or np.isnan(m):  #if there is no sowing month value
        result=[0]*51  
    else:                    #there is sowing moth value
        m=str(int(m)).zfill(2)
        if c=='' or np.isnan(c):            #no cultivar being identified
            result=group[[col for col in group.columns if 's'+m in col]].mean(axis=1).tolist()  #average across the three groups
        else:                              #month and cultivar are provided
            c=str(int(c)).replace("-",".")
            result=group.loc[:,'m'+c+'s'+m].tolist()
    return(result+[a])

def ProdBenefitsWithSegPara(crop,fer,G1):
    #production with reported month
    in_dir = RESULTS_ROOT / crop / fer
    yld=pd.read_csv(SIMOUT_DIR / f"EPIC_mz_result_20crv3_obs_{fer}.csv") 
    grid=pd.read_csv(INPUT_DIR / "African5minGrid_SIMUNIT.csv")[['SIMUNIT','POINT_X','POINT_Y']]
    grid.columns=['SIMUNIT','x','y'] 
    area=pd.read_csv(INPUT_DIR / "Africa_SIMUNIT_MZ_PhysicalArea.csv")  #area
    reportedsow=pd.read_csv(INPUT_DIR / "Africa_SimGrid_Confirmed_5min_4calibration.csv")[['SIMUNIT','ReportedSow_se1m','ReportedSow_se2m']]
    reportedsow=reportedsow.merge(area,how='left')
    reportedsow.columns=['SIMUNIT','s1m','s2m','PhysicalArea']
    #remove season 2 sowing month if it overlaps with the season 1, the difference of the two months is less than 4.
    reportedsow.loc[(reportedsow.s1m-reportedsow.s2m)%12<4,'s2m']=np.nan
    result=pd.DataFrame(columns=['season','nneigh','nseg','prod_mean','prod_cv'])
    for se in [1,2]:
        #se=1
        df = pd.read_csv(in_dir / f"{crop}_se{se}_{fer}_dis_withmon{G1}.csv")
        df.columns=['x','y','c','m']
        #only keep the simunits with simulated values
        temp1=df.merge(grid,how='inner',on=['x','y'])[['SIMUNIT','c','m']]
        temp1.columns=['SIMUNIT','sim_c','sim_m']
        temp1=temp1.groupby('SIMUNIT').agg(mode).reset_index()
        #First estimate the production with reported sow and mean over the three cultivars
        temp=reportedsow[['SIMUNIT','s'+str(se)+'m','PhysicalArea']]  #reported month and area
        temp.columns=['SIMUNIT','m','A']
        temp=temp.merge(temp1,how='inner').reset_index()
        temp.insert(1, 'c', '')                                #inset cultivar columns
        df1=yld.merge(temp,how='left')                          #merge to yld dataframe
        temp=df1.groupby('SIMUNIT').apply(IdentifyCulGroup,include_groups=False).reset_index()   #compute yearly yield
        df1=pd.DataFrame([[temp.iloc[i,1][j] for j in range(len(temp.iloc[i,1]))] for i in temp.index]).fillna(0)  #conver the dataframe
        prod=np.dot(df1.iloc[:,51].T,df1.iloc[:,:51])
        result.loc[len(result),]=[se,0,0]+[prod.mean(),prod.std()*100/prod.mean()]
        for nneigh in range(6,32,2):
            for nseg in range(1,100):
                df1=clustering(df,nneigh,nseg)
                #compute production for each segmentation strategy
                temp=pd.read_csv(INPUT_DIR / "African5minGrid_SIMUNIT.csv")[['SIMUNIT','POINT_X','POINT_Y']]
                temp.columns=['SIMUNIT','x','y'] 
                df1=df1.merge(temp,how='left')[['SIMUNIT','c','m','seg']]  #add SIMUNIT
                df1=df1.groupby(['SIMUNIT']).agg(mode).reset_index()
                temp=reportedsow[['SIMUNIT','PhysicalArea']].dropna()
                df1=df1.merge(temp,on='SIMUNIT',how='left')[['SIMUNIT','seg','c','m','PhysicalArea']]
                df1.columns=['SIMUNIT','seg','c','m','A']
                df2=yld.merge(df1,how='left')
                temp=df2.groupby('SIMUNIT').apply(IdentifyCulGroup,include_groups=False).reset_index()   #compute yearly yield
                df1=pd.DataFrame([[temp.iloc[i,1][j] for j in range(len(temp.iloc[i,1]))] for i in temp.index]).fillna(0)  #conver the dataframe
                prod=np.dot(df1.iloc[:,51].T,df1.iloc[:,:51])
                result.loc[len(result),]=[se,nneigh,nseg]+[prod.mean(),prod.std()*100/prod.mean()]
    result.to_csv(in_dir / "nmarket_temp_sea1.csv",index=False)



#Step 2 - Save segmentation as tif files - convert segmentation to three rasters with the value of month, cul-1.tiff, cul0.tiff, and cul1.tiff
def csv2tif(csvinput,outtif,lon_col='x',lat_col='y',value_col='value'):
    """
    Convert a csv file with latitude, longitude and value columns to a GeoTIFF
    """
    df=csvinput#pd.read_csv(csvinput)
    lons=df[lon_col].values
    lats=df[lat_col].values
    values=df[value_col].values
    df=df.dropna() #remove rows with NaN values in latitude, longitude or value
    min_lon,max_lon=np.min(lons),np.max(lons)
    min_lat,max_lat=np.min(lats),np.max(lats)
    resolution=0.1 #define resolution
    rows=int((max_lat-min_lat)/resolution)+1 #create the output raster
    cols=int((max_lon-min_lon)/resolution)+1
   
    raster=np.full((rows,cols),np.nan) #Initialize with NaN values
    #Populate the raster based on dataframe values
    for lat,lon,value in zip(lats,lons,values):
        col=int((lon-min_lon)/resolution)
        row=int((max_lat-lat)/resolution)
        if 0<=row<rows and 0<=col<cols: #check bounds
            raster[row,col]=value
    transform=from_origin(min_lon,max_lat,resolution,resolution) #create transform
    #Write the raster to a GeoTIFF file
    with rasterio.open(
            outtif,
            'w',
            driver='GTiff',
            height=raster.shape[0],
            width=raster.shape[1],
            count=1,
            dtype=rasterio.float32, #Use float32 to stor NaN values
            crs='EPSG:4326',     #Geographic coordinate system
            transform=transform
            ) as dst:
                dst.write(raster,1)

#season 1 neib=18, nmarket=25, season 2 neib=16, nmarket=5
def Seg2tiff(crop,fer,G1):
    para=[[18,25],[16,5]]  #34 market for season 1, and 8 market for season 2
    in_dir = RESULTS_ROOT / crop / fer
    for se in [1,2]:
        df = pd.read_csv(in_dir / f"{crop}_se{se}_{fer}_dis_withmon{G1}.csv") 
        df.columns=['x','y','c','m']
        # 假设列: lon, lat, variety, month
        # 2. 特征预处理
        # 提取属性特征进行聚类
        X_attr = df[['c', 'm']].values
        # 提取空间坐标
        coords = df[['x', 'y']].values
        # 标准化属性特征（防止月份的1-12比品种的-1-1权重过大）
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_attr)
        # 3. 创建空间连通性约束 (这是保证“空间相连”的关键)
        # 我们定义每个点只与其最近的 6 个邻居点连通（模拟蜂窝状或网格状布局）
        connectivity = kneighbors_graph(coords, n_neighbors=para[se-1][0], include_self=False)
        # 4. 执行带约束的层级聚类
        # n_clusters: 你想要分割出的市场个数
        # linkage: 'ward' 最小化方差，使每个区域内部属性最接近
        model = AgglomerativeClustering(
            n_clusters=para[se-1][1], 
            connectivity=connectivity, 
            linkage='ward'
        )
        df['seg'] = model.fit_predict(X_scaled)
        temp=df[['seg','c','m']].groupby(['seg']).agg(mode).reset_index()
        df=df[['x','y','seg']].merge(temp,how='left')
        for cul in df.c.unique():
            temp=df.loc[df['c']==cul,]
            csv2tif(temp,in_dir / f"{crop}_se{se}_cul{cul}_{fer}_dommon.tiff",lon_col='x',lat_col='y',value_col='m')

#Step 3 - Estimate Production Benefits for Countries
def EstimateCountryBenefits(crop,fer,G1):
    in_dir = RESULTS_ROOT / crop / fer
    para=[[18,25],[16,5]]  #One of the optimum segmentation strategy
    grid=pd.read_csv(INPUT_DIR / "African5minGrid_SIMUNIT.csv")[['SIMUNIT','POINT_X','POINT_Y']]
    grid.columns=['SIMUNIT','x','y']
    grid['x']=grid['x'].round(2)
    grid['y']=grid['y'].round(2)
    reportedsow=pd.read_csv(INPUT_DIR / "Africa_SimGrid_Confirmed_5min_4calibration.csv")[['SIMUNIT','ReportedSow_se1m','ReportedSow_se2m','country','A']]
    reportedsow.columns=['SIMUNIT','s1m','s2m','country','A']
    #remove season 2 sowing month if it overlaps with the season 1, the difference of the two months is less than 4.
    reportedsow.loc[(reportedsow.s1m-reportedsow.s2m)%12<4,'s2m']=np.nan
    yld=pd.read_csv(SIMOUT_DIR / f"EPIC_mz_result_20crv3_obs_{fer}.csv")
    result=pd.DataFrame(columns=['season','country','nmarket','w_wo_cluster']+['year_'+str(y) for y in range(1971,2022)]) #output results
    for se in [1,2,3]:  #2 is the second season but with only current reported potential area,3 is the second month with all potential areas
        df = pd.read_csv(in_dir / f"{crop}_se{se//2+1}_{fer}_dis_withmon{G1}.csv")
        df.columns=['x','y','c','m']
        #only keep the simunits with simulated values
        temp1=df.merge(grid,how='inner',on=['x','y'])[['SIMUNIT','c','m']]
        temp1.columns=['SIMUNIT','sim_c','sim_m']
        temp1=temp1.groupby('SIMUNIT').agg(mode).reset_index()
        #First estimate the production with reported sow and mean over the three cultivars
        temp=reportedsow[['SIMUNIT','s'+str(se//2+1)+'m','A','country']]  #reported month and area
        temp.columns=['SIMUNIT','m','A','country']
        temp=temp.merge(temp1,how='inner').reset_index()
        temp.insert(1, 'c', '')                                #inset cultivar columns
        df1=yld.merge(temp,how='left')                          #merge to yld dataframe
        for coun in df1.country.unique():
            temp=df[df1['country']==coun]
            if temp.shape[0]>0:
                temp=temp.groupby('SIMUNIT').apply(IdentifyCulGroup).reset_index()   #compute yearly yield
                temp=pd.DataFrame([[temp.iloc[i,1][j] for j in range(len(temp.iloc[i,1]))] for i in temp.index]).fillna(0)  #conver the dataframe
                temp=np.dot(temp.iloc[:,51].T,temp.iloc[:,:51])
                result.loc[len(result),]=[se,coun,0,'wocluster']+temp.tolist()
 
        #production with the clustering
        temp=clustering(df,para[se//2][0],para[se//2][1])  #cluster
        temp['x']=temp['x'].round(2)
        temp['y']=temp['y'].round(2)
        temp=temp.merge(grid,how='inner',on=['x','y'])[['SIMUNIT','c','m','seg']]
        df1=temp.groupby(['SIMUNIT']).agg(mode).reset_index() #.astype('int')  #get value by seg num
    
        #grids have not been clustered
        temp=reportedsow[~reportedsow['SIMUNIT'].isin(df1['SIMUNIT'].tolist())][['SIMUNIT','s'+str(se//2+1)+'m']]
        temp.columns=['SIMUNIT','m']
        temp.insert(1, 'c', '')
        temp.insert(1,'seg','')
        df=pd.concat([df1, temp], ignore_index=True)
        
        #merge area
        if se<3:
            temp=reportedsow[['SIMUNIT','s'+str(se%2+1)+'m','A','country']].dropna()  #using reported season area
        else:
            temp=reportedsow[['SIMUNIT','A','country']].dropna() #using all potential area
        df=df.merge(temp[['SIMUNIT','A','country']],how='left')
        df=yld.merge(df,how='left')
        for coun in df.country.unique():
            temp=df[df['country']==coun]  #find country
            seg=len(temp.seg.unique())
            if temp.shape[0]>0:
                temp=temp.groupby('SIMUNIT').apply(IdentifyCulGroup).reset_index()  #get yld
                temp=pd.DataFrame([[temp.iloc[i,1][j] for j in range(len(temp.iloc[i,1]))] for i in temp.index]).fillna(0)
                temp=np.dot(temp.iloc[:,51].T,temp.iloc[:,:51])  #estimate production
                result.loc[len(result),]=[se,coun,seg,'wcluster']+temp.tolist()
    result.to_csv(in_dir / f"{crop}_CountryProdBenefitByClustering_{fer}{G1}.csv", index=False)

#Step 4 - Plotting
def MainText_Plot_4(crop,fer,G1):
    #colors for the number of market
    cmap1 = plt.cm.tab20
    colors=cmap1.colors
    market=[25,5]
    africa=gpd.read_file(BOUNDARY_FILE)
    africa=africa[africa['CONTINENT']=='Africa']
    in_dir = RESULTS_ROOT / crop / fer
    fig = plt.figure(figsize=(18,15))  #,height_ratios=[2,1]
    gs=fig.add_gridspec(2,2,height_ratios=[2,1.2])
    #Benefits maps at the bottom
    color_cmap=['Greens','Purples','Oranges']
    nmarket=[] 
    for se in range(2):
        #Market distribution #process data
        ax_top=fig.add_subplot(gs[0,se])
        africa.boundary.plot(ax=ax_top,linewidth=0.5,color='black')
        #read data
        for c in [-1,0,1]:
            tif=in_dir / f"{crop}_se{se+1}_cul{c}_{fer}_dommon.tiff"
            if os.path.exists(tif):
                with rasterio.open(tif) as src:
                    data = src.read(1)  # 读取第一个波段数据
                    bounds=src.bounds  #get the bounds of the raster
                    #data_list.append(data)
                    nodata_value=src.nodata  #get the NoData value if it exists
                    ax_top.imshow(data,cmap=color_cmap[c+1],vmin=1,vmax=12,extent=[bounds.left,bounds.right,bounds.bottom,bounds.top])
        ax_top.set_xlim([-17.5,51.5])
        ax_top.set_ylim([-34.5,37.5])
        ax_top.set_xlabel("Longitude", fontsize=18)
        ax_top.set_ylabel("Latitude", fontsize=18)
        ax_top.set_xticks([-10,0,10,20,30,40,50],[r'$10^\circ$W',r'$0^\circ$',r'$10^\circ$E',r'$20^\circ$E',r'$30^\circ$',r'$40^\circ$E',r'$50^\circ$E'],fontsize=14)
        ax_top.set_yticks([-30,-20,-10,0,10,20,30],[r'$30^\circ$S',r'$20^\circ$S',r'$10^\circ$S',r'$0^\circ$',r'$10^\circ$N',r'$20^\circ$N',r'$30^\circ$N'],fontsize=14)
        ax_top.text(-28,41,chr(97+se),fontsize=26)
        ax_top.text(-16,-33, "Number of markets: "+ str(market[se])+"\n",style='italic',fontsize=15) 
    ########################################Country production benefits############################################
    ax_bottom=fig.add_subplot(gs[1,:])
    bar_width = 0.61
    bar_location=[0.1,0.965]
    tree_location=[[0.161,0.11,0.290,0.23],[0.574,0.11,0.290,0.172]]
    seg=pd.read_csv(in_dir / f"mz_CountryProdBenefitByClustering_{fer}{G1}.csv")
    for se in [1,3]:
        df=seg.groupby(['season','w_wo_cluster'])[seg.columns[4:]].sum().reset_index()
        x=df.loc[(df.season==se)&(df.w_wo_cluster=='wcluster'),df.columns[2:]].values-df.loc[(df.season==se)&(df.w_wo_cluster=='wocluster'),df.columns[2:]].values
        ax_bottom.bar(bar_location[se//2],(x/1000000).mean(),yerr=np.std(x/1000000),capsize=15,width=bar_width,color='white',edgecolor='black',label='Primary season')
        df=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wcluster')]
        df['prod_mean_wc']=df.iloc[:,4:].mean(axis=1)
        df=df[['country','nmarket','prod_mean_wc']]
        temp=seg.loc[(seg.season==se)&(seg.w_wo_cluster=='wocluster')]
        temp['prod_mean_woc']=temp.iloc[:,4:].mean(axis=1)
        temp=temp[['country','prod_mean_woc']]
        df=df.merge(temp,how='left')
        df['benefit']=(df['prod_mean_wc']-df['prod_mean_woc'])/1000000
        data=df[df['benefit']>=1].sort_values(by='benefit',ascending=False)
        temp=df[df['benefit']<1]
        data.loc[len(data)]=['Others',temp['nmarket'].mean(),temp['prod_mean_wc'].sum(),temp['prod_mean_woc'].sum(),temp['benefit'].sum()]
        data['benefit%']=(data['prod_mean_wc']-data['prod_mean_woc'])*100/data['prod_mean_woc']
        #change cell
        data.loc[data['country']=='United Republic of Tanzania','country']='Tanzania'
        data.loc[data['country']=='Democratic Republic of the Congo','country']='DRC'
        if se==1:
            data['benefit%']=data['benefit%'].round(0).astype(int)
            data.loc[data['benefit%']>200,'benefit%']=200
            data['benefit%']=data['benefit%'].astype('str')
            data.loc[data['benefit%']=="200",'benefit%']=">200"
        data['nmarket']=data['nmarket'].round(0).astype('int')
        nmarket=nmarket+data['nmarket'].tolist()
        se_ax=fig.add_axes(tree_location[se//2])
        labels=data['country']+"\n"+data['benefit'].round(1).astype('str')+"MT("+data['benefit%'].astype('str')+"%)"
        if se==3:
            labels=data['country']+"\n"+data['benefit'].round(1).astype('str')+"MT"
        squarify.plot(sizes=data['benefit'],alpha=0.7,label=labels,color=[colors[n] for n in data['nmarket']],
                      edgecolor='white',linewidth=1,text_kwargs={'fontsize':10},ax=se_ax)
        se_ax.axis('off')
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    maturity=['Early', 'Medium', 'Late']
    for i in range(3):
        cbar_ax = fig.add_axes([0.91, 0.444+0.145*i, 0.05, 0.145])
        cmap=mpl.cm.get_cmap(color_cmap[i],13)
        bounds=np.arange(13)
        norm=mpl.colors.BoundaryNorm(bounds,cmap.N)
        cb=mpl.colorbar.ColorbarBase(
            cbar_ax,cmap=cmap,norm=norm,boundaries=bounds,
            ticks=[],spacing='proportional',orientation='vertical'
        )
        cb.ax.tick_params(which='minor', color='white', labelcolor='white')
        cb.set_label(maturity[i],fontsize=14,rotation=270,labelpad=15)
        cbar_ax.yaxis.set_ticks([])
        #cbar_ax.set_frame_on(False)
        for j, label in enumerate(labels):
            y=0.03+(1/12)*j
            cbar_ax.text(0.35,y,label,va='center',ha='left',fontsize=9,transform=cbar_ax.transAxes)
    ax_bottom.text(-0.38,55,'c',fontsize=26)
    ax_bottom.set_ylabel("Production benefit (MT)", fontsize=18,labelpad=25)
    ax_bottom.set_xticks([0.1,0.95],['Primary season','Minor season'],fontsize=18)
    ax_bottom.yaxis.set_tick_params(labelsize=14)
    #colorbar for number of market
    nmarket=list(set(nmarket))
    cbar_ax = fig.add_axes([0.91,0.11,0.05,0.263])
    nmarket_cmap=mcolors.ListedColormap([colors[n] for n in nmarket])
    norm=mcolors.BoundaryNorm(boundaries=np.arange(0,len(nmarket)+1),ncolors=len(nmarket))
    sm=plt.cm.ScalarMappable(cmap=nmarket_cmap,norm=norm)
    sm.set_array([])
    cbar=fig.colorbar(sm,orientation="vertical",cax=cbar_ax)
    cbar.ax.tick_params(which='minor', color='white', labelcolor='white')
    cbar.set_ticks([0.5+x for x in range(len(nmarket))])
    cbar.set_ticklabels([f'{i}' for i in nmarket])
    cbar.set_label("Number of markets",fontsize=14,rotation=270,labelpad=25)
    fig.savefig(PLOTS_DIR / "MainText_Fig_4.png", format="png", dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "MainText_Fig_4.pdf", format="pdf", dpi=300, bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":
    crop = "mz"
    fer = "wfer_gridcalibratedWaHi"
    G1 = "_GT1"
    # ProdBenefitsWithSegPara(crop, fer, G1)  # Step 1
    # Seg2tiff(crop, fer, G1)  # Step 2
    # EstimateCountryBenefits(crop, fer, G1)  # Step 3
    MainText_Plot_4(crop, fer, G1)  # Step 4
