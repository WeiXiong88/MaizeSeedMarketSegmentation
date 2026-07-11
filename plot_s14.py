##########################################################FOR SUPPLEMENTARY FIGURE S14#######################################
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy import stats
except ImportError:
    stats = None


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
DATA_DIR = Path(r"D:\Works\AfricaMzSg")
INPUT_DIR = DATA_DIR / "data"
OUTPUT_DIR = DATA_DIR / "output"
RESULTS_ROOT = DATA_DIR / "results"


def linear_regression(x, y):
    if stats is not None:
        return stats.linregress(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    return slope, intercept, np.nan, np.nan, np.nan


def yld2prodInAtable(crop,fer,cropmodel):
    in_dir = INPUT_DIR
    result = pd.DataFrame(columns=['sowtype', 'type', 'cul', 'sea'] + [x for x in range(1971, 2022)])
    df = pd.read_csv(in_dir / "FAOSTAT_African_Countries_data_en_1-21-2025.csv")
    temp = df.loc[(df['Area'] == 'Africa') & (df['Year'] > 2018) & (df['Element'] == 'Production'), 'Value'].tolist()
    result.loc[len(result),] = ['fao', 'prod', 'real', '2019-2023'] + [np.nan] * (51 - len(temp)) + [x / 1000000 for x in temp]
    #########################################Estimated production by detrended FAO yield######################
    df = df[(df['Area'] == 'Africa') & (df['Year'] > 1970) & (df['Year'] < 2022)][['Element', 'Year', 'Value']]
    temp = df.loc[df['Element'] == 'Yield', ['Year', 'Value']]
    temp.columns = ['year', 'value']
    slope, intercept, r_value, p_value, std_err = linear_regression(temp['year'], temp['value'])
    yld = ((temp['value'] + (2023 - temp['year']) * slope) / 1000).tolist()  # detrend yield to 2020 level
    result.loc[len(result),] = ['fao', 'yld', 'yld_detrend2020', ''] + [round(x, 3) for x in yld]
    prod = (np.array(yld) * df[(df['Year'] > 2011) & (df['Year'] < 2022) & (df['Element'] == 'Area harvested')][
        'Value'].mean() / 1000000).tolist()
    result.loc[len(result),] = ['fao', 'prod', 'yid_detrend2020', ''] + prod
    ####################################SIMULATED YIELD#########################################################
    # Resport sowing month and area for the primary and minor maize reported sowing month
    sow = pd.read_csv(in_dir / "Africa_SimGrid_Confirmed_5min_4calibration.csv")[
        ['SIMUNIT', 'ReportedSow_se1m', 'ReportedSow_se2m', 'A']]
    sow.columns = ['SIMUNIT', 'se1m', 'se2m', 'A']
    area = pd.read_csv(in_dir / "Africa_SIMUNIT_MZ_PhysicalArea.csv")
    area = area.merge(sow[['SIMUNIT', 'A']], how='outer').fillna(0)
    area['se1a'] = area['PhysicalArea']  # Area of the primary maize is set to the physical area of maize in around 2020 (SPAM)
    area['se2a'] = area['A'] - area['PhysicalArea']  # Actual area of second maize is not know, setting to Physical-Harvest>0
    area[area < 0] = 0
    sow = sow[['SIMUNIT', 'se1m', 'se2m']]
    area = area[['SIMUNIT', 'se1a', 'se2a']]
    #######################################SIMULATED VALUE# WITHOUT FERTILIZER#################################
    sowing_window = ['ReportedSow', 'FixedSow', 'OptimumFixedSow', 'OptimumYearSow']
    in_dir = OUTPUT_DIR
    for sw in sowing_window:
        df = pd.read_csv(in_dir / (cropmodel+"_"+crop+"_yield_" + sw + "_" + fer + ".txt"), sep=r'\s+')  # change here for fertilizer
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

def Supplementary_Plot_S14(crop,fer,G1):
    cmap1 = plt.cm.tab20
    colors=cmap1.colors
    fig, ax = plt.subplots(2, 1, figsize=(12,12))
    in_dir = RESULTS_ROOT / crop / fer
    seg=pd.read_csv(in_dir / f"mz_CountryProdBenefitByClustering_{fer}{G1}.csv") 
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
    fig.savefig(PLOTS_DIR / "plot_s14.png",format="png",dpi=300,bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s14.pdf",format="pdf",dpi=300,bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    # Helper script: import yld2prodInAtable from here, or call it manually with crop, fertilizer, and model.
    crop='mz'
    fer="wfer_gridcalibratedWaHi"
    G1="_GT1"
    Supplementary_Plot_S14(crop,fer,G1)
