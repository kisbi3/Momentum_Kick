import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import multiprocessing
import time

import class_pool as classes

time_start = time.time()
'''numpy array print options'''
# np.set_printoptions(precision=5, suppress=True, linewidth = 205)
'''
    N_ch 에 대한 N_k 계산(using class)
'''

# 0부터가 아니라 -0.8 ~ 0.8로 해야될거 같은데??

#함수 선언



if __name__ == '__main__':
    # try:
    #     impact_param_initial = np.loadtxt("Nch.csv", usecols=[0])
    # except:
    #     impact_param_initial = None
    # Nch = np.arange(0, 131.5, 0.5)      # 0 ~ 131
    # multiple = 100000                     # fitting error가 너무 작아서 제대로 fitting이 안 되는것 처럼 보여 x100하여 에러 크게 하기

    # ''' 한번에 fitting '''
    # Nch_fitting = classes.FITTING(Nch, impact_param_initial, multiple)
    # impact_param, Nch = Nch_fitting.fitting()
    # time_fitting = time.time()
    # print(f"Fitting End : {time_fitting-time_start:.3f}")
    # '''impact param, Nch'''
    # np.savetxt("Nch.csv", list(zip(impact_param, Nch)))

    ''' 몇개씩 나눠서 fitting '''
    # result_Nch = np.array([])
    # result_b = np.array([])
    # length = 2
    # for i in range(int(len(impact_param_initial)/length)):
    #     impact_param_separate = impact_param_initial[length*i:length*i+length]
    #     Nch_separate = Nch[length*i : length*i+length]
    #     Nch_fitting = classes.FITTING(Nch_separate, impact_param_separate, multiple)
    #     impact_param, Nch_result = Nch_fitting.fitting()
    #     print("\n\nImpact parameter Result\n", impact_param)
    #     print("Multiplicity Result\n", Nch_result)
    #     result_Nch = np.append(result_Nch, Nch_result)
    #     result_b = np.append(result_b, impact_param)


    #     if i == int(len(impact_param_initial)/length)-1:
    #         impact_param_separate = impact_param_initial[length*i+length::]
    #         Nch_separate = Nch[length*i+length::]
    #         print(impact_param_separate)
    #         print(Nch_separate)
    #         Nch_fitting = classes.FITTING(Nch_separate, impact_param_separate, multiple)
    #         impact_param, Nch_result = Nch_fitting.fitting()
    #         print("\n\nImpact parameter Result\n", impact_param)
    #         print("Multiplicity Result\n", Nch_result)
    #         result_Nch = np.append(result_Nch, Nch_result)
    #         result_b = np.append(result_b, impact_param)
    
    '''아예 따로 fitting'''
    # Nch_error = 0
    # for i in range(len(impact_param_initial)):
    #     bb = impact_param_initial[i]
    #     multi = Nch[i]
    #     print("Initial: ", bb, multi)
    #     Nch_fitting = classes.FITTING(multi, bb, multiple)
    #     bb_result, multi_result = Nch_fitting.fitting()
    #     print("result: ", bb_result, multi_result)
    #     result_Nch = np.append(result_Nch, multi_result)
    #     result_b = np.append(result_b, bb_result)
    #     Nch_error = Nch_error + (multi_result - multi)*(multi_result - multi)
    


    # time_fitting = time.time()
    # print(f"Fitting End : {time_fitting-time_start:.3f}")
    # '''impact param, Nch'''
    # np.savetxt("Nch.csv", list(zip(result_b, result_Nch)))

    # print("\n\n Total Result")
    # print("\nImpact parameter : \n", result_b)
    # print("\nMultiplicity : \n", result_Nch)
    # print("Total Error: ", Nch_error)
    # # 4.80598032e-09 -> 3.92938671e-09 -> 7.95110233e-10 -> 7.00893618e-10 -> 6.97643659e-10
    # if impact_param_initial is None:
    #     pass
    # else:
    #     # print("\n Impact parameter fitting difference \n", impact_param - impact_param_initial)
    #     '''separate'''
    #     print("\n Impact parameter fitting difference \n", result_b - impact_param_initial)     
    ''' 꽤 정밀하게 완료. 이제 fitting 안할 예정(Error : 6.97643685e-10)'''
    
    ''' Calculate Nch '''
    # impact_parameter = np.loadtxt("Nch.csv", usecols=[0])
    # Nch = np.loadtxt("Nch.csv", usecols=[1])
    # calculate = classes.function_gpu()
    # Nk_calculate = calculate.Nk_b(impact_parameter)
    # print(Nk_calculate)

    # np.savetxt("bNchNk_result.csv", list(zip(impact_parameter, Nch, Nk_calculate)))

    ''' Drawing Graph '''
    Nch = np.loadtxt("bNchNk_result.csv", usecols=[1])
    Nk = np.loadtxt("bNchNk_result.csv", usecols=[2])
    
    fig = plt.figure()
    ax = plt.axes()
    fig.set_size_inches(35, 16.534, forward=True)    
    
    plt.plot(Nch, Nk, linewidth=7)
    # plt.legend()
    # plt.show()

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')
    plt.grid(color='silver',linestyle=':',linewidth=5)
    plt.xlabel(r'$ N_{ch} $',size=70)
    plt.ylabel(r'$ N_k $',size=70)
    plt.tight_layout()

    fig.savefig('Nch_Nk.png')


    time_end = time.time()
    print("time : ",time_end-time_start, "sec")