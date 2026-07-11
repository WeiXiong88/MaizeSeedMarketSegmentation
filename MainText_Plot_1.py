"""
Plot 1 in the mainTEXT - Simulated vs. reported Yield and production 
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

INPUT_DIR = Path(r"D:\Works\AfricaMzSg\input")
OUTPUT_DIR = Path(r"D:\Works\AfricaMzSg\output")

##########################################################FOR MAIN TEXT FIGURE 1#######################################
def yld2prodInAtable(fer):
    """
    A table containing area-weighted yield and production
    :param fer:
    :return: result
    """
    in_dir = INPUT_DIR
    # sea 1=primary 2=minor 3=total/combined
    result = pd.DataFrame(columns=['sowtype', 'type', 'cul', 'sea'] + [x for x in range(1971, 2022)])
    ###FAO Reported yield and production, computed from detrended yield, production = detrend yield  x max(area in 1971-2021)
    df = pd.read_csv(in_dir / "FAOSTAT_African_Countries_data_en_1-21-2025.csv")
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
        df = pd.read_csv(in_dir / f"EPIC_mz_yield_{sw}_{fer}.txt", sep=r'\s+')  # change here for fertilizer
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

def MainText_Plot_1(fer):
    """
    Variation of simulated vs. reported maize yield over year.
    """
    type_colors = {-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"}  # blue, red, green
    # result=pd.read_csv("D:\\works\\AfricaMzSg\\results\\mz\\"+fer+"\\1_mz_Africa_WeightedYld_Production_"+fer+G1+".csv")
    result = yld2prodInAtable(fer)
    type_colors = ["#00CD6C", "#AF58BA", "#FFC61E"]
    width_ratios = [1, 1]
    label = ['Early', 'Medium', 'Late']
    #####################################
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={'width_ratios': width_ratios})  #
    sw = ['fao', 'ReportedSow', 'FixedSow','OptimumYearSow']  # ['ReportedSow','FixedSow','OptimumFixedSow','OptimumDecadeSow','OptimumYearSow']
    label_name = ['Simulated Early', 'Simulated Medium', 'Simulated Late']
    x = list(range(1971, 2022))
    y = result[(result['sowtype'] == 'fao') & (result['type'] == 'yld')].iloc[0, 4:].values
    axes[0].plot(x, y, color='grey', label="Reported (FAO)")  # FAO
    for cul in range(3):
        y1 = result[(result['sowtype'] == 'ReportedSow') & (result['type'] == 'yld') & (result['sea'] == 3) & (
                    result['cul'] == (cul + 1))].iloc[0, 4:].values * 0.7
        y2 = result[(result['sowtype'] == 'OptimumYearSow') & (result['type'] == 'yld') & (result['sea'] == 3) & (
                    result['cul'] == (cul + 1))].iloc[0, 4:].values * 0.7
        y = result[(result['sowtype'] == 'OptimumFixedSow') & (result['type'] == 'yld') & (result['sea'] == 3) & (
                    result['cul'] == (cul + 1))].iloc[0, 4:].values * 0.7
        axes[0].fill_between(x, [float(y) for y in y1], [float(y) for y in y2], color=type_colors[cul], alpha=0.1)
        axes[0].plot(x, y, color=type_colors[cul], label=label_name[cul])
    # Bar
    # mzarea=area[['se1a','se2a']].sum().sum()
    xtick = [1, 2, 3, 4]
    bar_width = 0.66
    # fao
    # mean=np.mean(result[(result['sowtype']=="fao")&(result['type']=="yld")][[str(yr) for yr in range(1971,2022)]].values*mzarea/1000000)
    # sd=np.std(result[(result['sowtype']=="fao")&(result['type']=="yld")][[str(yr) for yr in range(1971,2022)]].values*mzarea/1000000)
    mean = result.loc[(result['sowtype'] == "fao") & (result['cul'] == "real"), [yr for yr in range(2017, 2022)]].mean(axis=1).tolist()[0]
    sd = result.loc[(result['sowtype'] == "fao") & (result['cul'] == "real"), [yr for yr in range(2017, 2022)]].std(axis=1).tolist()[0]
    axes[1].bar(xtick[0], mean, bar_width, yerr=sd, color='grey', error_kw=dict(lw=1, capsize=2, capthick=1))
    # simulated
    for ty in [1, 2, 3]:  # sowing type
        for cul in [1, 2, 3]:  # variety
            mean1 = result[(result['sowtype'] == sw[ty]) & (result['type'] == "prod") & (result['sea'] == 1) & (result['cul'] == cul)][
                [yr for yr in range(1971, 2022)]].values  # *(area.se1a.sum())/1000000
            mean2 = result[(result['sowtype'] == sw[ty]) & (result['type'] == "prod") & (result['sea'] == 2) & (result['cul'] == cul)][
                [yr for yr in range(1971, 2022)]].values  # *(area.se2a.sum())/1000000
            sd = (mean1 + mean2).std()
            axes[1].bar(xtick[ty] - bar_width / 3 + (bar_width / 3) * (cul - 1), mean2.mean(), bar_width / 3, color=type_colors[cul - 1])
            axes[1].bar(xtick[ty] - bar_width / 3 + (bar_width / 3) * (cul - 1), mean1.mean(), bar_width / 3,
                        color=type_colors[cul - 1], yerr=sd,
                        error_kw=dict(lw=1, capsize=2, capthick=1))
    reportedsow_mean = result.loc[(result.sowtype == "ReportedSow") & (result.type == "prod"), ['cul'] + [yr for yr in
                                                                                                          range(1971,
                                                                                                                2022)]].groupby(
        'cul').sum().mean().mean()
    # 104
    axes[1].bar(2, reportedsow_mean, bar_width, edgecolor='black', linewidth=1, facecolor='none')
    axes[0].set_xlabel("Year", fontsize=15)
    axes[0].set_yticks([1.5, 2, 2.5, 3, 3.5])
    axes[0].set_ylabel("Yield (t/ha)", fontsize=15)
    axes[1].set_ylabel("Production (MT)", fontsize=15)
    axes[1].set_xticks([1, 2, 3, 4],
                       ['Reported\nFAO', 'Simulated with\nreported sowing', 'Simulated with\nidentical\noptimum sowing',
                        'Simulated with\nvaried optimum\n sowing'], fontsize=10)
    axes[0].text(1962, axes[0].get_ylim()[1] + (axes[0].get_ylim()[1] - axes[0].get_ylim()[0]) / 20, "a",
                 fontsize=20)  # 3.82
    axes[1].text(-0.10, axes[1].get_ylim()[1] + (axes[1].get_ylim()[1] - axes[1].get_ylim()[0]) / 20, "b",
                 fontsize=20)  # 195
    axes[0].legend(loc='lower center', ncol=4, fontsize=10, bbox_to_anchor=(1.1, -0.30))
    axes[1].annotate(
        "Production \nestimated\nfrom FAO data",  # The text to display
        xy=(1, 98),  # The point the arrow points to
        xytext=(0.595, 150),  # The point where the text is placed
        arrowprops=dict(facecolor="black", arrowstyle="->", linewidth=0.5),  # Arrow properties
        fontsize=10
    )
    axes[1].text(1.5, 170, "Productions estimated from simulations\n with varied variety maturity", fontsize=10)
    plt.subplots_adjust(left=0.06, bottom=0.22, top=0.88, right=0.99)
    fig.savefig(PLOTS_DIR / "MainText_Fig_1.png", format="png", dpi=300)
    fig.savefig(PLOTS_DIR / "MainText_Fig_1.pdf", format="pdf", dpi=300)


if __name__ == "__main__":
    fer = "wfer_gridcalibratedWaHi"
    MainText_Plot_1(fer)  # Plot 1 in Main Text.
