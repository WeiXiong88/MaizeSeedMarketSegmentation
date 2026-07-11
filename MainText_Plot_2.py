"""
Plot 2 in the mainTEXT - Changes in SI 
"""
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
#%matplotlib inline

SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"D:\Works\AfricaMzSg\input")  #SCRIPT_DIR.parent / "data"
OUTPUT_DIR = Path(r"D:\Works\AfricaMzSg\output") #SCRIPT_DIR.parent / "output"

##########################################################FOR MAIN TEXT FIGURE 2#######################################
def AggregatedYld(fer):
    in_dir = OUTPUT_DIR
    sw=["ReportedSow","OptimumFixedSow"]
    #Maize physical area
    area=pd.read_csv(INPUT_DIR / "Africa_SimGrid_Confirmed_5min_4calibration.csv")[['SIMUNIT','ReportedSow_se1m','ReportedSow_se2m','A']]
    area.columns=['SIMUNIT','se1m','se2m','A']
    yld=pd.DataFrame(columns=['sowtype','matu','sea','GT1']+[str(yr) for yr in range(1971,2022)])
    for sow in sw:
        df0=pd.read_csv(in_dir / f"EPIC_mz_yield_{sow}_{fer}.txt",sep=r'\s+')  #Change here when you want to use wofer yield
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
#yld.to_csv("D:\\works\\AfricaMzSg\\results\\mz\\"+fer+"\\mz_aggregatedAfricanYield_"+fer+".csv", index=False)

def MainText_Plot_2(fer,season):
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
    ax1.set_ylabel('Optimum Sow            Reported Sow',fontsize=14)
    ax1.set_xlabel('Suitability Index (SI)',fontsize=14)
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
        ax2.set_ylabel("SI (t/ha)",fontsize=14,labelpad=10)
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
        ax3.set_xlabel("Decade",fontsize=14)
        ax3.set_ylabel("SI (t/ha)",fontsize=14)
        ax3.tick_params(axis='both', which='major')
    ax2.legend(bbox_to_anchor=(1.2, 1.03), loc='upper right',fontsize=10)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.30, hspace=0.15)
    ax1.text(-0.6,5.8,'a',fontsize=20)
    ax2.text(1968,ax2.get_ylim()[1]+(ax2.get_ylim()[1]-ax2.get_ylim()[0])/20,'b',fontsize=20)
    ax3.text(1968,ax3.get_ylim()[1]+(ax3.get_ylim()[1]-ax3.get_ylim()[0])/20,'c',fontsize=20)
    fig.savefig(PLOTS_DIR / "MainText_Fig_2.png",format="png",dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(PLOTS_DIR / "MainText_Fig_2.pdf",format="pdf",dpi=300, bbox_inches='tight', pad_inches=0)

if __name__ == "__main__":
    fer = "wfer_gridcalibratedWaHi"
    MainText_Plot_2(fer, season=1)  # Plot 2 in Main Text.
