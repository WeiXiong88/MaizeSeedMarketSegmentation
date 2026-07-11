from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import ListedColormap


TRAIL_DIR = "../#"##Path(__file__).resolve().parent
TRIAL_SITE_FILE = "../Data/MaizeMaturity_Phenotype_Mean_2011-2016_4_Simlation.csv"
fig_dir = "../Plots/" #str(TRAIL_DIR) + "/"


def load_trial_sites():
    trial = pd.read_csv(TRIAL_SITE_FILE)
    trial = trial.dropna(subset=["Latitude", "Longitude", "Year"]).copy()
    trial_sites = (
        trial.groupby(["Latitude", "Longitude"], as_index=False)
        .agg(
            experimental_years=("Year", "nunique"),
            n_records=("Year", "size"),
            country=("Country", lambda x: "; ".join(sorted(set(map(str, x))))),
            mega_environment=("MegaEnviron", lambda x: "; ".join(sorted(set(map(str, x))))),
        )
    )
    return trial_sites


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
    trial_sites = load_trial_sites()
    marker_sizes = 18 + trial_sites["experimental_years"] * 18
    ax.scatter(
        trial_sites["Longitude"],
        trial_sites["Latitude"],
        s=marker_sizes + 28,
        facecolors="none",
        edgecolors="white",
        linewidths=2.8,
        zorder=5,
    )
    ax.scatter(
        trial_sites["Longitude"],
        trial_sites["Latitude"],
        s=marker_sizes,
        facecolors="none",
        edgecolors="#0072B2",
        linewidths=1.4,
        zorder=5,
        label="Trial sites",
    )
    # cbar = fig.colorbar(img, ax=ax, ticks=np.arange(1,13),orientation='vertical', fraction=0.1, pad=0.05)
    patches = [
        plt.plot([], [], marker="s", ms=10, ls="", color=colors[i], label=class_name[i + 1])[0]
        for i in range(len(colors))
    ]
    year_levels = [1, 3, int(trial_sites["experimental_years"].max())]
    year_levels = sorted(set(year_levels))
    trial_handles = [
        plt.scatter(
            [],
            [],
            s=18 + y * 18,
            facecolors="none",
            edgecolors="#0072B2",
            linewidths=1.4,
            label=f"{y} experimental year{'s' if y > 1 else ''}",
        )
        for y in year_levels
    ]
    legend1 = ax.legend(
        handles=patches,
        title="Maize Mega-Environments",
        bbox_to_anchor=(1.01, 1.01),
        loc='upper left',
        frameon=True
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=trial_handles,
        title="Trial sites",
        bbox_to_anchor=(1.01, 0.43),
        loc="upper left",
        frameon=True,
        scatterpoints=1,
    )
    # "viridis"
    # ax.set_xlim([-17.62504269,51.13387])
    # ax.set_ylim(lat_range)
    ax.set_xlabel("Longitude", fontsize=18)
    ax.set_ylabel("Latitude", fontsize=18)
    ax.set_xticks(
        [-20, -10, 0, 10, 20, 30, 40, 50],
        [r'$20^\circ$W', r'$10^\circ$W', r'$0^\circ$', r'$10^\circ$E',
         r'$20^\circ$E', r'$30^\circ$', r'$40^\circ$E', r'$50^\circ$E'],
        fontsize=14,
    )
    ax.set_yticks(
        [-30, -20, -10, 0, 10, 20, 30, 40],
        [r'$30^\circ$S', r'$20^\circ$S', r'$10^\circ$S', r'$0^\circ$',
         r'$10^\circ$N', r'$20^\circ$N', r'$30^\circ$N', r'$40^\circ$N'],
        fontsize=14,
    )

    #plt.show()
    fig.savefig(fig_dir + "plot_s2.png", format="png", dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(fig_dir + "plot_s2.pdf", format="pdf", dpi=300, bbox_inches='tight', pad_inches=0)
Supplementary_Plot_S2()
