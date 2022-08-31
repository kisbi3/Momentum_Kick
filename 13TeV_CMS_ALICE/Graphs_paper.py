import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit


mpl.rcParams["text.usetex"] = True

# alice 논문에서는 near-side 정의를 -1.28<phi<1.28으로 두었다.
# fig.set_size_inches(50, 30, forward=True)
ali_dat=[]
ali_phi=[]
ali_err=[]
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=129, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=167, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=205, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=7))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=129, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=167, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=205, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=7))
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=129, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=129, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=129, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=129, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=167, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=167, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=167, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=167, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=205, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=205, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=205, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=205, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=7)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=7)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[6],skiprows=12,max_rows=7)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[7],skiprows=12,max_rows=7)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)

#pt 1~4를 해야하나?
cms_dat=[]
cms_phi=[]
cms_err=[]
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[1],skiprows=14, max_rows= 9))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[0],skiprows=14, max_rows= 9))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=9)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=9)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=9)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=9)
cms_err.append((err_sta1**2+err_sys1**2)**0.5)
cms_err.append((err_sta2**2+err_sys2**2)**0.5)

atl_dat=[]
atl_phi=[]
atl_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_atlas/atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[1],skiprows=2,max_rows=14))
atl_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_atlas/atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[0],skiprows=2,max_rows=14))

#ALICE CMS ALICE CMS ... 순서
mom_dat=[]
mom_phi=[]
mom_dat.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt0-5.csv',delimiter=',',usecols=[3],skiprows=1))
mom_dat.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[2],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt0-5.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[0],skiprows=1))

# 적분이 조금 특수해서 직접 계산한 후에 넣어줘야 할듯
# ali_integral = np.trapz(ali_dat[3], x = ali_phi[3])
ali_integral = 0.016963325
# cms_integral = np.trapz(cms_dat[3], x = cms_phi[3])
# cms_integral = 0.026520954
cms_integral = 0.0257493594
# mom_ali_integral = np.trapz(mom_dat[6], x = mom_phi[3])
# mom_cms_integral = np.trapz(mom_dat[7], x = mom_phi[3])
mom_ali_integral = np.loadtxt('pTdis_integral.csv', delimiter=',', usecols = [1])
mom_cms_integral = np.loadtxt('pTdis_integral.csv', delimiter=',', usecols = [3])
cms_ratio = cms_integral/mom_cms_integral
ali_ratio = ali_integral/mom_ali_integral
# print(ali_integral, cms_integral, mom_ali_integral, mom_cms_integral, ali_ratio, cms_ratio)

fig1, axes1 = plt.subplots(nrows=1, ncols=4,figsize=(100,20))

#그래프 그리기
for i in range(4):
        if i == 3:
                axes1[i].plot(mom_phi[i], mom_dat[2*i] - min(mom_dat[2*i]), color = "blue", linewidth=7, linestyle='-')
                axes1[i].errorbar(atl_phi[0], atl_dat[0]-min(atl_dat[0]), color="blue", markersize=30, marker='o', linestyle=' ', linewidth=5, capsize=15)
                axes1[i].set_title(r'$0.5<p_{T,\,trig(assoc)}<5$', size = 70, pad=30)
                axes1[i].set_xlabel(r'$\Delta\phi$', size=70)
        
        # elif i==4:
        #         axes1[i].plot(mom_phi[i], mom_dat[2*i-1]*ali_ratio, color = "red", linewidth=7, linestyle='-')
        #         axes1[i].plot(mom_phi[i], mom_dat[2*i]*cms_ratio, color = "black", linewidth=7, linestyle='-')
        #         axes1[i].errorbar(ali_phi[i-1], ali_dat[i-1], yerr=(abs(ali_err[2*(i-1)+1]),ali_err[2*(i-1)]), color="red", markersize=20, marker='o', linestyle=' ', fillstyle='none', linewidth=5, capsize=15)
        #         axes1[i].errorbar(cms_phi[i-1], cms_dat[i-1], yerr=(abs(cms_err[2*(i-1)+1]),cms_err[2*(i-1)]), color="black", markersize=20, marker='o', linestyle=' ', fillstyle='none', linewidth=5, capsize=15)
        #         axes1[i].set_xlabel(r'$p_{T,\,trig(assoc)}$', size=70)
        #         axes1[i].set_ylabel(r'$Y^{ridge}$', size=70)
        else: 
                axes1[i].plot(mom_phi[i], mom_dat[2*i] - min(mom_dat[2*i]), color = "red", linewidth=7, linestyle='-')
                axes1[i].plot(mom_phi[i], mom_dat[2*i+1] - min(mom_dat[2*i+1]), color = "black", linewidth=7, linestyle='-')
                axes1[i].errorbar(ali_phi[i], ali_dat[i]-min(ali_dat[i]), yerr=(abs(ali_err[2*i+1]),ali_err[2*i]), color="red", linestyle=' ', linewidth=7, capthick=3, capsize=15)
                axes1[i].errorbar(cms_phi[i], cms_dat[i]-min(cms_dat[i]), yerr=(abs(cms_err[2*i+1]),cms_err[2*i]), color="black", linestyle=' ', linewidth=7, capthick=3, capsize=15)
                axes1[i].scatter(ali_phi[i], ali_dat[i]-min(ali_dat[i]), edgecolors="red", s=600, marker='o', facecolors='none', linewidths=7)
                axes1[i].scatter(cms_phi[i], cms_dat[i]-min(cms_dat[i]), edgecolors="black", s=600, marker='o', facecolors='none', linewidths=7)
                st = i+1
                en = i+2
                axes1[i].set_title(str(st)+r'$<p_{T,\,trig(assoc)}<$'+str(en), size = 70, pad=30)
                # axes1[i].set_ylim(-0.001,0.018)
                # axes1[i].tick_params(labelleft = False)
                axes1[i].set_xlabel(r'$\Delta\phi$', size=70)
                # axes1[i].tight_layout()
        axes1[i].minorticks_on()
        # axes1[0].tick_params(labelleft = True)

        # axes1[i].set_ylim(-0.001,0.02)
        # axes1[i].set_xlim(-1.28,1.28)

        axes1[i].tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right='true')
        axes1[i].tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right='true')
        axes1[i].grid(color='silver',linestyle=':',linewidth=3)

# axes1[i].legend(framealpha=False, fontsize = 70)
axes1[0].set_ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$', size=70)
fig1.tight_layout(h_pad=-1)

fig1.savefig('./paper_graph/Paper_phiCorr.png')
fig1.savefig('./paper_graph/Paper_phiCorr.pdf')

fig2, axis2 = plt.subplots(nrows=1, ncols=1,figsize=(40,20))

i=4
# print(mom_dat[2*i-1])
# print(mom_dat[2*i])
# correction_pt = np.loadtxt('pTdis_correction.csv',delimiter=',',usecols=[0],skiprows=1)
# ali_correction = np.loadtxt('pTdis_correction.csv',delimiter=',',usecols=[1],skiprows=1)
# cms_correction = np.loadtxt('pTdis_correction.csv',delimiter=',',usecols=[2],skiprows=1)
# corr_index = np.where(mom_phi[i] == correction_pt)
# ali_ratio = ali_correction/mom_dat[2*i-1][corr_index]
# cms_ratio = cms_correction/mom_dat[2*i][corr_index]
# print(mom_phi[i])
# print(mom_dat[2*i-1]*ali_ratio)
# print(ali_ratio)
# print()
axis2.plot(mom_phi[i], mom_dat[2*i-1]*ali_ratio, color = "red", linewidth=7, linestyle='-')
axis2.plot(mom_phi[i], mom_dat[2*i]*cms_ratio, color = "black", linewidth=7, linestyle='-')
# axis2.scatter(mom_phi[i], mom_dat[2*i-1]*ali_ratio, color = "red")
# axis2.scatter(mom_phi[i], mom_dat[2*i]*cms_ratio, color = "black")
axis2.errorbar(ali_phi[i-1], ali_dat[i-1], yerr=(abs(ali_err[2*(i-1)+1]),ali_err[2*(i-1)]), color="red", linestyle=' ', linewidth=7, capthick=3, capsize=15)
axis2.errorbar(cms_phi[i-1], cms_dat[i-1], yerr=(abs(cms_err[2*(i-1)+1]),cms_err[2*(i-1)]), color="black", linestyle=' ', linewidth=7, capthick=3, capsize=15)
axis2.scatter(ali_phi[i-1], ali_dat[i-1], edgecolors="red", s=600, marker='o', facecolors='none', linewidths=7)
axis2.scatter(cms_phi[i-1], cms_dat[i-1], edgecolors="black", s=600, marker='o', facecolors='none', linewidths=7)
axis2.set_xlabel(r'$p_{T,\,trig(assoc)}$', size=70)
axis2.set_ylabel(r'$Y^{ridge}$', size=70)
axis2.minorticks_on()
axis2.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right='true')
axis2.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right='true')
axis2.grid(color='silver',linestyle=':',linewidth=3)
fig2.tight_layout()
fig2.savefig('./paper_graph/Paper_pTdis.png')
fig2.savefig('./paper_graph/Paper_pTdis.pdf')


fig3 = plt.figure()
ax = plt.axes()
fig3.set_size_inches(35, 16.534, forward=True)

x = np.arange(0,4,0.01)

f1 = 0.32*np.exp(0.93*x)        #7TeV frnk
# f2 = 8*np.exp(-3/x)        #13TeV frnk
f2 = 0.83+0.5*x*x
f3 = np.zeros(len(x))           #200GeV
f22 = 0.66*np.exp(0.71*x)
PbPb_fr = np.exp(-1.395/x)
PbPb_nk = 20.2*np.exp(-0.207*x)
f4 = PbPb_fr*PbPb_nk            #PbPb 2.76TeV frnk

# plt.plot(x,f1, color = 'red', linewidth=7, label=r'$pp, \, 7TeV$')
plt.plot(x,f3+4, color = 'red', linewidth=7, label=r'$AuAu, \, 200GeV$')
# plt.plot(x,f3, color = 'red', linewidth=7, label=r'$0.66e^{0.71p_T}$')
plt.plot(x,f4, color = 'black', linewidth=7, label=r'$PbPb, \, 2.76TeV$')
plt.plot(x,f2, color = 'blue', linewidth=7, label=r'$pp, \, 13TeV$')
plt.plot(x,f22, color = 'blue', linewidth=7, linestyle='--', label=r'$pp, \, 13TeV$')
# plt.plot(x,f2, color = 'blue', linewidth=7, label=r'$0.83+0.5{p_T}^2$')

plt.xlabel(r'$p_T^{trig}$',size=70)
plt.ylabel(r'$f_{R} \langle N_k \rangle $',size=70)
# plt.yscale('log')

plt.xlim(0,4)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')

plt.grid(color='silver',linestyle=':',linewidth=5)
plt.legend(fontsize=45, loc='upper left')

plt.tight_layout()

fig3.savefig('./paper_graph/Paper_frnk.png')
fig3.savefig('./paper_graph/Paper_frnk.pdf')


# fig3.clear()
# plt.plot(x, 2*np.exp(-2/x), label='1')
# plt.plot(x, 2*np.exp(-1/x), label='2')
# plt.plot(x, 1*np.exp(-2/x), label='3')
# plt.legend()
# plt.show()