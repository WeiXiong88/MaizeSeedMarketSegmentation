##########################################################FOR SUPPLEMENTARY FIGURE S7#######################################
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy import stats
except ImportError:
    stats = None


SCRIPT_DIR = Path(__file__).resolve().parent
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"
INPUT_DIR = Path(r"D:\works\AfricaMzSg\input")
OUTPUT_DIR = Path(r"D:\works\AfricaMzSg\output")


def linear_regression(x, y):
    if stats is not None:
        return stats.linregress(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    return slope, intercept, np.nan, np.nan, np.nan


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
    slope, intercept, r_value, p_value, std_err = linear_regression(temp['year'], temp['value'])
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
        df = pd.read_csv(in_dir / ("mz_yield_" + sw + "_" + fer + ".txt"), sep=r"\s+")  # change here for fertilizer
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
def Supplementary_Plot_S8(fer):
    type_colors = {-1: "#00CD6C", 0: "#AF58BA", 1: "#FFC61E"}  # blue, red, green
    result = yld2prodInAtable(fer)
    result = result[(result.sowtype != 'fao') & (result.type == 'yld') & (result.sea != 3)]
    result['mean'] = result.iloc[:, 4:].mean(axis=1)
    type = ['Planting', 'Season', 'Maturity', 'Year']
    data = []  # Create a empty list
    data.append(result[['sowtype', 'mean']].groupby('sowtype').mean().reset_index()['mean'].tolist())  # Planting
    data.append(result[['sea', 'mean']].groupby('sea').mean().reset_index()['mean'].tolist())  # Season
    data.append(result[['cul', 'mean']].groupby('cul').mean().reset_index()['mean'].tolist())  # Maturity
    data.append(result[[yr for yr in range(1971, 2022)]].mean().tolist())  # Year
    fig = plt.figure(figsize=(6, 4.5))
    plt.boxplot(data, vert=True, patch_artist=True, showmeans=False, showfliers=False,
                boxprops=dict(facecolor='lightblue', color="blue"),
                medianprops=dict(color="red", linewidth=2),
                whiskerprops=dict(color="blue", linewidth=1.5),
                capprops=dict(color="blue", linewidth=1.5))
    plt.xticks([1, 2, 3, 4], type, fontsize=12)
    plt.ylabel("Mean yield (t/ha)", fontsize=15)
    #plt.show()
    fig.savefig(PLOTS_DIR / "plot_s8.png",format="png",dpi=300, bbox_inches='tight')
    fig.savefig(PLOTS_DIR / "plot_s8.pdf",format="pdf",dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    Supplementary_Plot_S8("wfer_gridcalibratedWaHi")
