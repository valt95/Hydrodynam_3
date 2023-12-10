# Функции для расчета излучения
import numpy as np
# Различные физические константы
k = 1.38e-23 # Постоянная Больцмана Дж/К
m_p = 1.67e-27 # Масса протона кг
c = 3e8 # Скорость света
pk_to_meter = 3e16 # Переводной коэффициент пк->м
G = 6.67e-11 # гравитационная постоянная Н*м^2/кг^2
M_sun = 2e30 # масса солнца в кг
gamma = 5/3 # коэффициент адиабаты
def eddington_light(n_H, x, E, M_pbh, r_phys, M_unit):  # Не учитывается рекомбинация так как светит во все стороны (соответственно зависимости от температуры здесь нет)
    # Количество свободных электронов
    # x = np.ones(x.size)
    n_e = x * n_H
    n_H_0 = (1-x) * n_H
    # Количество нейтрального водорода

    L_edd = 10**40 * M_pbh/1e9

    E_0 = 13.6   # в эВ
    sigma_photoionization = 6e-22   # м**2
    k_ion = n_H_0 * sigma_photoionization * (E/E_0)**(-3) * (E>E_0) #(с обрезкой меньше 13.6 эВ)
    sigma_scat = 6.652e-29  # м**2

    H = np.zeros(r_phys.size)
    H[0] = L_edd/(4*np.pi*(r_phys[0]*pk_to_meter)**2)
    Pressure_euler_coupling = np.zeros(r_phys.size-1)
    Energy_euler_coupling = np.zeros(r_phys.size-1)
    for i in range(1,r_phys.size-1):
        H[i] = (r_phys[i-1]/r_phys[i])**2 * H[i-1] * (1 - (r_phys[i]-r_phys[i-1]) * pk_to_meter * (k_ion[i-1] + n_e[i-1]*sigma_scat))
        if H[i]<0:
            H[i] = 0
        elif H[i]>H[0]:
            H[i] = 0
        elif H[i]>H[i-1]:
            H[i] = 0
        Pressure_euler_coupling[i-1] = -H[i-1] * (k_ion[i-1] + n_e[i-1]*sigma_scat)
        Energy_euler_coupling[i-1] = -H[i-1] * (k_ion[i-1] + n_e[i-1]*sigma_scat)
    # H = H*np.where(H>0,1,0)
    # H = H*np.where(H<=H[0],1,0)
    # for i in range(1,r_phys.size-1):
    #     if H[i]>H[i-1]:
    #         H[i] = 0

    # Energy_euler_coupling = (H[1:] - H[:-1])*4*np.pi
    # Pressure_euler_coupling = (H[1:] - H[:-1])*(4*np.pi/c)
    Energy_euler_coupling *= 4*np.pi*np.where(H[:-1]!=0,1,0)
    Pressure_euler_coupling *= 4*np.pi/c*np.where(H[:-1]!=0,1,0)
    return H, Pressure_euler_coupling/M_unit, Energy_euler_coupling/M_unit

def ioniz_frac_easy(dt, x, n_H, E, T, H):
    # Темп ионизации
    E_0 = 13.6
    n_H_0 = (1-x)*n_H
    n_e = x*n_H
    sigma_photoionizatioh = 6e-22 # м**2
    ioniz_rate = H * n_H_0*sigma_photoionizatioh*(E/E_0)**(-3)*(E>E_0)  # 1/сек*м**3

    # Темп рекомбинации
    recombination_coef = 2.6*10**(-13)*(T/10**4)**(-0.85) * 1e-6 # м**3/cек значение из статьи Кара
    recombination_rate = n_e**2 * recombination_coef  # 1/сек*м**3

    # Вычисляем приросты
    n_H_0 += dt*(recombination_rate - ioniz_rate)
    n_e += dt*(ioniz_rate - recombination_rate)

    # print('recomb',recombination_rate[0])
    # print('ioniz',ioniz_rate[0])
    for i in range(n_e.size):
        if n_e[i]<0:
            n_e[i] = 0
            n_H_0[i] = n_H[i]
        elif n_e[i]>n_H[i]:
            n_e[i] = n_H[i]
            n_H_0[i] = 0

        if n_H_0[i]<0:
            n_H_0[i] = 0
            n_e[i] = n_H[i]
        elif n_H_0[i]>n_H[i]:
            n_H[i] = n_H_0
            n_e[i] = 0

    x = n_e/n_H

    return x

def ioniz_frac(H, n_H, x, E, r_phys, T, M_pbh, M_unit,dt):
    # Темп ионизации
    E_0 = 13.6
    sigma_photoionizatioh = 6e-22  # м**2
    alpha = 2.6*10**(-13)*(T/10**4)**(-0.85) * 1e-6 # м**3/c
    ioniz_rate = np.zeros(x.size)
    recomb_rate = np.zeros(x.size)
    time_of_light_sim = 0
    for j in range(x.size-1):
        dt_light = (r_phys[j+1]-r_phys[j])*3e16/c # время за которое свет проходит данную ячейку
        ioniz_rate[j] = sigma_photoionizatioh * (E_0/E)**3 * H[j] * n_H[j] * (1-x[j])
        recomb_rate[j] = alpha[j] * n_H[j]**2 * x[j]**2
        x_actual = np.zeros(2)     # Переменная для расчета степени ионизации
        x_actual[0] = 1
        x_actual[1] = x[j]
        while np.abs(x_actual[1]-x_actual[0])>0.01*x_actual[1]:  # Реализация численного решения уравнения для степени ионизации
            x_actual[0] = x_actual[1]
            x_actual[1] = x_actual[0] + (ioniz_rate[j] - recomb_rate[j])/(n_H[j])*dt_light
            # if j==0:
            #     print('0', x_actual[0], 'j=', j)
            #     print('1', x_actual[1])
            #     print('ioniz',ioniz_rate[j])
            #     print('recomb', recomb_rate[j])
            # print('multiplier',dt_light/(n_H[j]))
            # print(dt_light*3e-8)
            # Блок условий который удерживает x в рамках [0,1]
            if x_actual[1]<0:
                x_actual[1] = 0
                break
            elif x_actual[1]>1:
                x_actual[1] = 1
                break
            time_of_light_sim += dt_light
            # Условие, чтобы вычисления не выходили за пределы основного шага по времени
            if time_of_light_sim>=dt:
                break
        x[j] = x_actual[1] # Новое равновесное значение x в ячейке
        # print(x[j])
        # Пересчет потока излучения (функция не оптимизирована, считаются два ненужных здесь значения + лишние пространственные ячейки)
        H, dump_1, dump_2 = eddington_light(n_H, x, E, M_pbh, r_phys, M_unit)



    return x

def easy_eddington_light(x, r_phys, M_pbh):
    L_edd = 10**40 * M_pbh/(1e9)
    sigma_scat = 6.652e-29  # м**2
    H = L_edd/(4*np.pi*(r_phys*pk_to_meter)**2)
    easy_Pressure_euler_coupling = x[:-1]*H[:-1]*sigma_scat/(c*m_p)
    easy_Energy_euler_coupling = 0
    return H, easy_Pressure_euler_coupling, easy_Energy_euler_coupling