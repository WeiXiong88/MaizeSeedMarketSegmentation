from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"D:\works\AfricaMzSg\input")
BOUNDARY_FILE = Path(
    r"E:\RSG Dropbox\Wei Xiong\Works\CurrentProcessing\0_AfricanMaizeSorghum\Plots"
    r"\CountryBoundaryLines\ne_110m_admin_0_countries.shp"
)


def Supplementary_Plot_S7():
    africa=gpd.read_file(BOUNDARY_FILE)
    africa=africa[africa['CONTINENT']=='Africa']
    #Plot for comparison between reported and simulated maize yields
    df0=pd.read_csv(INPUT_DIR / "mz_calibrated.csv")
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
                           [r'$20^\circ$W',r'$10^\circ$W',r'$0^\circ$',r'$10^\circ$E',r'$20^\circ$E',r'$30^\circ$',
                            r'$40^\circ$E',r'$50^\circ$E'],fontsize=14)
        axes[r,0].set_yticks([-30,-20,-10,0,10,20,30,40],
                           [r'$30^\circ$S',r'$20^\circ$S',r'$10^\circ$S',r'$0^\circ$',r'$10^\circ$N',r'$20^\circ$N',
                            r'$30^\circ$N',r'$40^\circ$N'],fontsize=14)
        #sc1 = axes[0].scatter(df['lon'],df['lat'], c=df['obsy'], cmap='viridis')
        #africa.boundary.plot(ax=axes[0],linewidth=0.5,color='black')
        sc2 = axes[r,1].imshow(raster2, extent=[lon_min, lon_max, lat_min, lat_max], origin='lower', cmap='viridis', alpha=0.8)
        africa.boundary.plot(ax=axes[r,1],linewidth=0.5,color='black')
        #fig.colorbar(sc2,ax=axes[r,1],label='Simulated Yields')
        axes[r,1].set_xlabel("Longitude", fontsize=18)
        axes[r,1].set_ylabel("Latitude", fontsize=18)
        axes[r,1].set_xticks([-20,-10,0,10,20,30,40,50],
                           [r'$20^\circ$W',r'$10^\circ$W',r'$0^\circ$',r'$10^\circ$E',r'$20^\circ$E',r'$30^\circ$',
                            r'$40^\circ$E',r'$50^\circ$E'],fontsize=14)
        axes[r,1].set_yticks([-30,-20,-10,0,10,20,30,40],
                           [r'$30^\circ$S',r'$20^\circ$S',r'$10^\circ$S',r'$0^\circ$',r'$10^\circ$N',r'$20^\circ$N',
                            r'$30^\circ$N',r'$40^\circ$N'],fontsize=14)
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
    fig.savefig(PLOTS_DIR / "plot_s7.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s7.pdf",format="pdf",dpi=300, bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    Supplementary_Plot_S7()
