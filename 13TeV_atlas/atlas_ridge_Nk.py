import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

mpl.rcParams["text.usetex"] = True

#atlas data들
phi_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[0])
data_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[1])
phi_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[0])
data_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[1])
phi_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[0])
data_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[1])
phi_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[0])
data_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[1])
phi_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[0])
data_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[1])
phi_13TeV_80_90=np.loadtxt('./atlasgraphs/13TeV_80~90.csv',delimiter=',',usecols=[0])
data_13TeV_80_90=np.loadtxt('./atlasgraphs/13TeV_80~90.csv',delimiter=',',usecols=[1])
phi_13TeV_70_80=np.loadtxt('./atlasgraphs/13TeV_70~80.csv',delimiter=',',usecols=[0])
data_13TeV_70_80=np.loadtxt('./atlasgraphs/13TeV_70~80.csv',delimiter=',',usecols=[1])
phi_13TeV_60_70=np.loadtxt('./atlasgraphs/13TeV_60~70.csv',delimiter=',',usecols=[0])
data_13TeV_60_70=np.loadtxt('./atlasgraphs/13TeV_60~70.csv',delimiter=',',usecols=[1])
phi_13TeV_50_60=np.loadtxt('./atlasgraphs/13TeV_50~60.csv',delimiter=',',usecols=[0])
data_13TeV_50_60=np.loadtxt('./atlasgraphs/13TeV_50~60.csv',delimiter=',',usecols=[1])
phi_13TeV_40_50=np.loadtxt('./atlasgraphs/13TeV_40~50.csv',delimiter=',',usecols=[0])
data_13TeV_40_50=np.loadtxt('./atlasgraphs/13TeV_40~50.csv',delimiter=',',usecols=[1])
phi_13TeV_30_40=np.loadtxt('./atlasgraphs/13TeV_30~40.csv',delimiter=',',usecols=[0])
data_13TeV_30_40=np.loadtxt('./atlasgraphs/13TeV_30~40.csv',delimiter=',',usecols=[1])
#계산값들
Nk_059_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.590000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_059_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.590000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_060_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.600000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_060_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.600000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_061_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.610000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_061_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.610000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_062_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.620000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_062_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.620000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_063_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.630000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_063_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.630000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_064_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.640000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_064_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.640000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_065_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.650000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_065_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.650000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_066_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.660000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_066_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.660000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_067_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.670000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_067_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.670000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_068_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.680000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_068_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.680000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_090_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.900000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_090_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.900000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_080_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.800000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_080_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk0.800000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_110_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk1.100000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_110_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk1.100000.csv',delimiter=',',usecols=[3], skiprows=1)
Nk_120_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk1.200000.csv',delimiter=',',usecols=[0], skiprows=1)
Nk_120_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_Nk1.200000.csv',delimiter=',',usecols=[3], skiprows=1)
# Multi_095_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[0], skiprows=1)
# Multi_095_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[3], skiprows=1)
# Multi_105_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[0], skiprows=1)
# Multi_105_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[3], skiprows=1)
# Multi_115_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[0], skiprows=1)
# Multi_115_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[3], skiprows=1)
# Multi_125_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[0], skiprows=1)
# Multi_125_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[3], skiprows=1)
# Multi_135_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[0], skiprows=1)
# Multi_135_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[3], skiprows=1)


Multi_035_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_35.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_035_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_35.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_045_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_45.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_045_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_45.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_055_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_55.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_055_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_55.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_065_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_65.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_065_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_65.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_075_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_75.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_075_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_75.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_085_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_85.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_085_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_85.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_095_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_095_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_105_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_105_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_115_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_115_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_125_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_125_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[3], skiprows=1)
Multi_135_phi=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_140.csv',delimiter=',',usecols=[0], skiprows=1)
Multi_135_atl=np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-5_140.csv',delimiter=',',usecols=[3], skiprows=1)


fig1, axes1 = plt.subplots(nrows=3, ncols=3,figsize=(90,90),sharey='row', sharex='col')


# axes1[0][0].scatter(phi_13TeV_40_50, data_13TeV_40_50-min(data_13TeV_40_50), color = 'blue', s=1000, marker='o')
axes1[0][0].scatter(phi_13TeV_50_60, data_13TeV_50_60-min(data_13TeV_50_60), color = 'blue', s=2000, marker='o')
axes1[0][1].scatter(phi_13TeV_60_70, data_13TeV_60_70-min(data_13TeV_60_70), color = 'blue', s=2000, marker='o')
axes1[0][2].scatter(phi_13TeV_70_80, data_13TeV_70_80-min(data_13TeV_70_80), color = 'blue', s=2000, marker='o')
axes1[1][0].scatter(phi_13TeV_80_90, data_13TeV_80_90-min(data_13TeV_80_90), color = 'blue', s=2000, marker='o')
axes1[1][1].scatter(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), color = 'blue', s=2000, marker='o')
axes1[1][2].scatter(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), color = 'blue', s=2000, marker='o')
axes1[2][0].scatter(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), color = 'blue', s=2000, marker='o')
axes1[2][1].scatter(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), color = 'blue', s=2000, marker='o')
axes1[2][2].scatter(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), color = 'blue', s=2000, marker='o')
#양쪽 끝 색은 파란색, 파란색으로 그래프 사이 채우기, 그리고 중간의 그래프는 빨강!
color_q = ['blue', 'blue']

# axes1[0][0].plot(Multi_045_phi, Multi_045_atl-min(Multi_045_atl), color = color_q[1], linewidth=7, linestyle = '-')
axes1[0][0].plot(Multi_055_phi, Multi_055_atl-min(Multi_055_atl), color = color_q[1], linewidth=14, linestyle = '-')
axes1[0][1].plot(Multi_065_phi, Multi_065_atl-min(Multi_065_atl), color = color_q[1], linewidth=14, linestyle = '-')
axes1[0][2].plot(Multi_075_phi, Multi_075_atl-min(Multi_075_atl), color = color_q[1], linewidth=14, linestyle = '-')
axes1[1][0].plot(Multi_085_phi, Multi_085_atl-min(Multi_085_atl), color = color_q[1], linewidth=14, linestyle = '-')

#90<N<100
axes1[1][1].plot(Multi_095_phi, Multi_095_atl-min(Multi_095_atl), color = color_q[1], linewidth=14, linestyle = '-')
# axes1[0].plot(Nk_059_phi, Nk_059_atl-min(Nk_059_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,<N_k> \,=\,0.59$')
# axes1[0].plot(Nk_060_phi, Nk_060_atl-min(Nk_060_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$<N_k> \,=\,0.60$')
# axes1[0].plot(Nk_061_phi, Nk_061_atl-min(Nk_061_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,<N_k> \,=\,0.61$')
# axes1[0].fill_between(Nk_059_phi, Nk_059_atl-min(Nk_059_atl), Nk_061_atl-min(Nk_061_atl), color=color_q[0], alpha=0.4)
#100<N<110
axes1[1][2].plot(Multi_105_phi, Multi_105_atl-min(Multi_105_atl), color = color_q[1], linewidth=14, linestyle = '-')
# axes1[1].plot(Nk_060_phi, Nk_060_atl-min(Nk_060_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,<N_k> \,=\,0.60$')
# axes1[1].plot(Nk_062_phi, Nk_062_atl-min(Nk_062_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$<N_k> \,=\,0.62$')
# axes1[1].plot(Nk_062_phi, Nk_062_atl-min(Nk_062_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,<N_k> \,=\,0.62$')
# axes1[0].fill_between(Nk_060_phi, Nk_060_atl-min(Nk_060_atl), Nk_062_atl-min(Nk_062_atl), color=color_q[0], alpha=0.4)
#110<N<120
axes1[2][0].plot(Multi_115_phi, Multi_115_atl-min(Multi_115_atl), color = color_q[1], linewidth=14, linestyle = '-')
# axes1[2].plot(Nk_061_phi, Nk_061_atl-min(Nk_061_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,<N_k> \,=\,0.61$')
# axes1[2].plot(Nk_090_phi, Nk_090_atl-min(Nk_090_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$<N_k> \,=\,0.90$')
# axes1[2].plot(Nk_063_phi, Nk_063_atl-min(Nk_063_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,<N_k> \,=\,0.63$')
# axes1[2].fill_between(Nk_061_phi, Nk_061_atl-min(Nk_061_atl), Nk_063_atl-min(Nk_063_atl), color=color_q[0], alpha=0.4)
#120<N<130
axes1[2][1].plot(Multi_125_phi, Multi_125_atl-min(Multi_125_atl), color = color_q[1], linewidth=14, linestyle = '-')
# axes1[3].plot(Nk_062_phi, Nk_062_atl-min(Nk_062_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,<N_k> \,=\,0.62$')
# axes1[3].plot(Nk_110_phi, Nk_110_atl-min(Nk_110_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$<N_k> \,=\,1.10$')
# axes1[3].plot(Nk_064_phi, Nk_064_atl-min(Nk_064_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,<N_k> \,=\,0.64$')
# axes1[3].fill_between(Nk_062_phi, Nk_062_atl-min(Nk_062_atl), Nk_064_atl-min(Nk_064_atl), color=color_q[0], alpha=0.4)
#130<N
axes1[2][2].plot(Multi_135_phi, Multi_135_atl-min(Multi_135_atl), color = color_q[1], linewidth=14, linestyle = '-')
# axes1[4].plot(Nk_063_phi, Nk_063_atl-min(Nk_063_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,<N_k> \,=\,0.63$')
# axes1[4].plot(Nk_120_phi, Nk_120_atl-min(Nk_120_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$<N_k> \,=\,1.20$')
# axes1[4].plot(Nk_065_phi, Nk_065_atl-min(Nk_065_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,<N_k> \,=\,0.65$')
# axes1[4].fill_between(Nk_063_phi, Nk_063_atl-min(Nk_063_atl), Nk_065_atl-min(Nk_065_atl), color=color_q[0], alpha=0.4)

# axes1[0][0].text(-0.9, 0.0185, r'$40 \leq N^{rec}_{ch} < 50$', size = 70)
axes1[0][0].text(-0.9, 0.0145, r'$50 \leq N^{rec}_{ch} < 60$', size = 150)
axes1[0][1].text(-0.9, 0.0145, r'$60 \leq N^{rec}_{ch} < 70$', size = 150)
axes1[0][2].text(-0.9, 0.0145, r'$70 \leq N^{rec}_{ch} < 80$', size = 150)
axes1[1][0].text(-0.9, 0.032, r'$80 \leq N^{rec}_{ch} < 90$', size = 150)
axes1[1][1].text(-0.9, 0.032, r'$90 \leq N^{rec}_{ch} < 100$', size = 150)
axes1[1][2].text(-0.9, 0.032, r'$100 \leq N^{rec}_{ch} < 110$', size = 150)
axes1[2][0].text(-0.9, 0.0625, r'$110 \leq N^{rec}_{ch} < 120$', size = 150)
axes1[2][1].text(-0.9, 0.0625, r'$120 \leq N^{rec}_{ch} < 130$', size = 150)
axes1[2][2].text(-0.9, 0.0625, r'$130 \leq N^{rec}_{ch}$', size = 150)

# axes1[0][0].text(-1.1, 0.008, r'$N_{trk}^{offline}<35$', size = 80)

axes1[0][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 150)
axes1[1][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 150)
axes1[2][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 150)
axes1[0][0].set_ylim(-0.001,0.0159)
axes1[1][0].set_ylim(-0.001,0.0349)
axes1[2][0].set_ylim(-0.001,0.0699)
for j in range(3):
    for i in range(3):
        if j==2:
            axes1[j][i].set_xlabel(r'$\Delta\phi$', size=150)

        axes1[j][i].minorticks_on()

        # axes1[j][i].set_ylim(-0.001,0.07)
        axes1[j][i].set_xlim(-1.1,1.1)

        axes1[j][i].tick_params(axis='both',which='major',direction='in',width=4,length=35,labelsize=100, top = 'true', right='true')
        axes1[j][i].tick_params(axis='both',which='minor',direction='in',width=4,length=20,labelsize=100, top = 'true', right='true')
        axes1[j][i].grid(color='silver',linestyle=':',linewidth=3)

        # axes1[i][j].legend(framealpha=False, fontsize = 70)

fig1.tight_layout(h_pad = -1)
fig1.savefig('./result/rezero_atlas_Nk.png')
# fig1.savefig('/home/jaesung/Desktop/Dropbox/논문/rezero_atlas_Nk.png')
# , transparent = True

fig1.clear()




# 일단 모든 데이터를 불러와보자.

cms_dat=[]
cms_phi=[]
cms_err1=[]
cms_err2=[]
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table1.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table3.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table5.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table7.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table9.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table11.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table13.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table15.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table17.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table19.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table21.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table23.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
# cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table25.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))

# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table1.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table3.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table5.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table7.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table9.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table11.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table13.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table15.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table17.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table19.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table21.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table23.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
# cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table25.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))

# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table1.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table1.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table3.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table3.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table5.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table5.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table7.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table7.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table9.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table9.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table11.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table11.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table13.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table13.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table15.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table15.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table17.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table17.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table19.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table19.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table21.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table21.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table23.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table23.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
# cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table25.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
# cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table25.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err1.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err2.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))


mom_dat=[]
mom_phi=[]
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_17.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_17.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_17.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_17.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_57.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_57.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_57.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_57.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_92.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_92.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_92.csv', delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_92.csv', delimiter=',',usecols=[2],skiprows=1))
# mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_155.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_155.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_155.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_155.csv',delimiter=',',usecols=[2],skiprows=1))

# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_17.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_17.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_17.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_17.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_57.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_57.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_57.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_57.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_92.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_92.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_92.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_92.csv',delimiter=',',usecols=[0],skiprows=1))
# mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt0-1_155.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt1-2_155.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt2-3_155.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('./phiCorr_Nk/phiCorrelation_pt3-4_155.csv',delimiter=',',usecols=[0],skiprows=1))


fig1, axes1 = plt.subplots(nrows=3, ncols=3,figsize=(105,105), sharey='row', sharex='col')
axes1[0][0].set_ylim(-0.0009,0.004)
axes1[1][0].set_ylim(-0.0015,0.01)
axes1[2][0].set_ylim(-0.001,0.018)
#그래프 그리기
for i in range(3):
    for j in range(3):
        # axes1[i][j].set_title(r'$0.5<p_{T,\,trig(assoc)}<5$', size = 70, pad=30)
        axes1[i][j].plot(mom_phi[3*i+j], mom_dat[3*i+j] - min(mom_dat[3*i+j]), color = "black", linewidth=7, linestyle='-')
        # axes1[i][j].plot(mom_phi[i], mom_dat[2*i+1] - min(mom_dat[2*i+1]), color = "black", linewidth=7, linestyle='-')
        axes1[i][j].errorbar(cms_phi[3*i+j], cms_dat[3*i+j]-min(cms_dat[3*i+j]), yerr=(abs(cms_err2[3*i+j]),cms_err1[3*i+j]), color="black", linestyle=' ', linewidth=14, capthick=6, capsize=30)
        axes1[i][j].scatter(cms_phi[3*i+j], cms_dat[3*i+j]-min(cms_dat[3*i+j]), edgecolors="black", s=2000, marker='o', facecolors='none', linewidths=14)
        # print(3*i+j, cms_dat[3*i+j], cms_dat[3*i+j]-min(cms_dat[3*i+j]))
        # print(3*i+j, cms_err1[3*i+j], cms_err2[3*i+j])

        
        st = j+1
        en = j+2
        if i==0:
            axes1[i][j].set_title(str(st)+r'$<p_{T,\,trig(assoc)}<$'+str(en), size = 200, pad=40)
        if i==2:
            axes1[i][j].set_xlabel(r'$\Delta\phi$', size=120)
        # axes1[i][j].tight_layout()
        axes1[i][j].minorticks_on()

        # axes1[i][j].set_ylim(-0.001,0.02)
        # axes1[i][j].set_xlim(-1.28,1.28)

        axes1[i][j].tick_params(axis='both',which='major',direction='in',width=4,length=35,labelsize=100, top = 'true', right='true')
        axes1[i][j].tick_params(axis='both',which='minor',direction='in',width=4,length=20,labelsize=100, top = 'true', right='true')
        axes1[i][j].grid(color='silver',linestyle=':',linewidth=3)


# axes1[0][0].text(-1.1, 0.008, r'$N_{trk}^{offline}<35$', size = 120)
axes1[0][0].text(-1.1, 0.0035, r'$35 \leq N_{trk}^{offline}<80$', size = 150)
axes1[1][0].text(-1.1, 0.009, r'$80 \leq N_{trk}^{offline}<105$', size = 150)
axes1[2][0].text(-1.1, 0.016, r'$105 \leq N_{trk}^{offline}$', size = 150)
# axes1[i][j].legend(framealpha=False, fontsize = 70)
axes1[0][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size=150)
axes1[1][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size=150)
axes1[2][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size=150)
# axes1[3][0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size=120)
fig1.tight_layout(h_pad=-1)

fig1.savefig('./result/cms_multiplicity_13TeV.png')

fig1.clear()


mom_Multix = np.loadtxt('./phiCorr_Nk/Multiplciity_dependence.csv', delimiter=',',usecols=[0],skiprows=1)
mom_Multiy = np.loadtxt('./phiCorr_Nk/Multiplciity_dependence.csv', delimiter=',',usecols=[1],skiprows=1)
mom_Multix_atl = np.loadtxt('./phiCorr_Nk/Multiplciity_dependence_atlas.csv', delimiter=',',usecols=[0],skiprows=1)
mom_Multiy_atl = np.loadtxt('./phiCorr_Nk/Multiplciity_dependence_atlas.csv', delimiter=',',usecols=[1],skiprows=1)

cms_Multi_x = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[0],skiprows=14)
cms_Multi_y = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[1],skiprows=14)
cms_Multi_sta1 = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[2],skiprows=14)
cms_Multi_sta2 = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[3],skiprows=14)
cms_Multi_sys1 = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[4],skiprows=14)
cms_Multi_sys2 = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table35.csv',delimiter=',',usecols=[5],skiprows=14)
multierr1 = (cms_Multi_sta1**2+cms_Multi_sys1**2)**0.5
multierr2 = (cms_Multi_sta2**2+cms_Multi_sys2**2)**0.5

atl_Multi_x=np.loadtxt('./atlasgraphs/Multiplicity_yield.csv',delimiter=',',usecols=[0], skiprows=1)
atl_Multi_y=np.loadtxt('./atlasgraphs/Multiplicity_yield.csv',delimiter=',',usecols=[1], skiprows=1)
# print(cms_Multi_x, cms_Multi_y)

fig1, axes1 = plt.subplots(figsize=(20,20))

axes1.text(2, 0.055, r'$1.0<p_T<2.0$ (GeV/c)', size = 50)
axes1.text(2, 0.052, r'$2.0<|\Delta\eta|<4.0$', size = 50)
axes1.text(2, 0.049, r'$0.5<p_T<5.0$ (GeV/c)', size = 50, color='blue')
axes1.text(2, 0.046, r'$2.0<|\Delta\eta|<5.0$', size = 50, color='blue')
axes1.set_ylabel('Associated Yield (GeV/c)', size = 70)
axes1.set_xlabel(r'$N_{trk}^{offline}, \,\, N_{ch}^{rec}$', size=70)

axes1.minorticks_on()

axes1.scatter(atl_Multi_x, atl_Multi_y, edgecolors="blue", s=600, marker='o', facecolors='none', linewidths=7)
axes1.errorbar(cms_Multi_x, cms_Multi_y, yerr=(multierr2, multierr1), color="black", linestyle=' ', linewidth=7, capthick=3, capsize=15)
axes1.scatter(cms_Multi_x, cms_Multi_y, edgecolors="black", s=600, marker='o', facecolors='none', linewidths=7)
axes1.plot(mom_Multix, mom_Multiy-min(mom_Multiy), color = "black", linewidth=7, linestyle='-')
axes1.plot(mom_Multix_atl, mom_Multiy_atl-min(mom_Multiy_atl), color = "blue", linewidth=7, linestyle='-')

axes1.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
axes1.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
axes1.grid(color='silver',linestyle=':',linewidth=3)


fig1.tight_layout(h_pad = -12)
fig1.savefig('./result/CMS_Multi_yield.png')
# fig1.savefig('/home/jaesung/Desktop/Dropbox/논문/rezero_atlas_Nk.png')
# , transparent = True

fig1.clear()
