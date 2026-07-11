from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats



SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"D:\works\AfricaMzSg\input")
OUTPUT_DIR = Path(r"D:\works\AfricaMzSg\output")

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
    slope, intercept, r_value, p_value, std_err = stats.linregress(temp['year'], temp['value'])
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
        inputfile = cropmodel+"_"+crop+"_yield_" + sw + "_" + fer + ".txt"
        df = pd.read_csv(in_dir / inputfile, sep=r'\s+')  # change here for fertilizer
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

def Supplementary_Plot_S15(crop,fer):
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
    fig.savefig(PLOTS_DIR / "plot_s15.png",format="png",dpi=300,bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "plot_s15.pdf",format="pdf",dpi=300,bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    Supplementary_Plot_S15("mz", "wfer_gridcalibratedWaHi")
