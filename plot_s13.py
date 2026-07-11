##########################################################FOR SUPPLEMENTARY FIGURE S12#######################################
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
RESULTS_ROOT = Path(r"D:\Works\AfricaMzSg\results")


def Supplementary_Plot_S13(crop,fer,G1):
    cmap1 = plt.cm.tab20
    colors=cmap1.colors
    fig, ax = plt.subplots(2, 1, figsize=(12,12))
    in_dir=RESULTS_ROOT / crop / fer
    seg=pd.read_csv(in_dir / ("mz_CountryProdBenefitByClustering_"+fer+G1+".csv")) 
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
    fig.savefig(PLOTS_DIR / "plot_s13.png",format="png",dpi=300,bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s13.pdf",format="pdf",dpi=300,bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    Supplementary_Plot_S13("mz", "wfer_gridcalibratedWaHi", "_GT1")
