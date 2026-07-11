##########################################################FOR SUPPLEMENTARY FIGURE S1#######################################
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap, Normalize


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"D:\works\AfricaMzSg\input")
TIF_DIR = INPUT_DIR / "05nc"
BOUNDARY_FILE = Path(
    r"E:\RSG Dropbox\Wei Xiong\Works\CurrentProcessing\0_AfricanMaizeSorghum\Plots"
    r"\CountryBoundaryLines\ne_110m_admin_0_countries.shp"
)


def Supplementary_Plot_S1():
    #Define color
    #type_colors = [{-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"},{-1: "red", 0: "lightgrey", 1: "green"}] #blue, red, green, maturity type
    #change_colors = {-1: "red", 0: "lightgrey", 1: "blue"} #change
    #African country boundary
    africa=gpd.read_file(BOUNDARY_FILE)
    africa=africa[africa['CONTINENT']=='Africa']
    tfile=["sow_month","har_area"]
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
        tif=TIF_DIR / ("MZ_season"+str(r%2+1)+"_"+tfile[r//2]+".tif")  #tif file for plotting
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
                           [r'$20^\circ$W',r'$10^\circ$W',r'$0^\circ$',r'$10^\circ$E',r'$20^\circ$E',r'$30^\circ$',
                            r'$40^\circ$E',r'$50^\circ$E'],fontsize=14)
        axes[r//2,r%2].set_yticks([-30,-20,-10,0,10,20,30,40],
                           [r'$30^\circ$S',r'$20^\circ$S',r'$10^\circ$S',r'$0^\circ$',r'$10^\circ$N',r'$20^\circ$N',
                            r'$30^\circ$N',r'$40^\circ$N'],fontsize=14)
        axes[r//2,r%2].text(-35,45,chr(97+r),fontsize=26)
    fig.savefig(PLOTS_DIR / "plot_s1.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s1.pdf",format="pdf",dpi=300, bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    Supplementary_Plot_S1()  # Reported sowing month and harvest area
