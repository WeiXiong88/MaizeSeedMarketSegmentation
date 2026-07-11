##########################################################FOR SUPPLEMENTARY FIGURE S9#######################################
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
RESULTS_ROOT = Path(r"D:\Works\AfricaMzSg\results")
para = [[18, 25], [16, 5]]


def Supplementary_Plot_S10(crop,fer,G1):
    point=[[151215226.6,2.708938],[48290267.12,5.104918]]  #season 1 neib=24, nmarket=31, season 2 neib=10, nmarket=4
    in_dir=RESULTS_ROOT / crop / fer
    df=pd.read_csv(in_dir / (crop+"_ProdBenefitsToSegmentationPara_"+fer+G1+".csv"))
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
    fig.savefig(PLOTS_DIR / "plot_s10.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s10.pdf",format="pdf",dpi=300, bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    Supplementary_Plot_S10("mz", "wfer_gridcalibratedWaHi", "_GT1")
