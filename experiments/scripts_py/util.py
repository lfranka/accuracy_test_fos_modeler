import rsf.api as rsf
import matplotlib.pyplot as plt
import numpy as np
import os

def loadrsf(FileName : str) -> np.ndarray:
    '''
        Exporta os valores de um arquivo .rsf para um numpy ndarray

        Entrada:
            - FileName -> Nome do arquivo .rsf
        Saindo:
            - data -> Vetor numpy com o tanho da primeiro dimensão do arquivo .rsf    
    '''

    if not FileName.endswith('.rsf'):
        FileName += '.rsf'

    inp     = rsf.Input(FileName)
    data    = np.zeros(inp.int('n1'), dtype='f')
    inp.read(data)

    return data


def convolve(a : np.ndarray, b : np.ndarray, nt : int, dt : float, ot : int) -> np.ndarray: 
    '''
        Faz a convolução entre dois dados e janela ele no tempo

        Entradas:
            a  -> Vetor para convolução
            b  -> Segundo vetor para convolução
            nt -> Número de amostras no tempo
            dt -> Espaçamento das amostras no tempo
            ot -> Origem das amostras no tempo
        
        Saida:
            data[mask] -> Vetor de convolução janelado com base
                          no eixo do tempo.
    '''
    data = np.convolve(a, b, mode='full')

    t_data = np.arange(len(data)) * dt

    mask = (t_data >= ot) & (t_data < (nt * dt))

    return data[mask]


def create_sap(G : np.ndarray, wav : np.ndarray, t : np.ndarray, 
               rho : float, dx : float, dz : float) -> np.ndarray:

    dwav = np.gradient(wav, t)

    nt = len(t)
    dt = t[1] - t[0]
    ot = t[0]

    sap = convolve(G, dwav, nt, dt, ot)

    fator = (-rho) * dt * dx * dz

    sap *= fator

    return sap


def diff_data(a : np.ndarray, b : np.ndarray) -> np.ndarray:
    
    dif  = ((a - b) / max(b)) * 100

    return dif


def load_data_cg(theta : np.ndarray, FileName : str):

    Dir = os.path.join(os.getcwd(), 'dados')

    # carregar dados teóricos
    cg_file = []
    for i in range(0, len(theta)):
        cg_tem = np.loadtxt(os.path.join(Dir, f'{FileName}_{i}.txt'), dtype='f')
        cg_file.append(cg_tem)

    cg_file = np.array(cg_file)

    f       = cg_file[0, :, 0]        # eixo da frequência teórica

    cg = []
    for i in range(0, len(theta)):
        cg_tem  = cg_file[i, :, 1]        # velocidade de grupo teórica
        cg.append(cg_tem)

    cg = np.array(cg)

    return f, cg


def calc_cg(theta : np.ndarray, f : np.ndarray, dt : float, dh : float, c : float) -> np.ndarray:

    pi = np.pi

    cgf = np.zeros((len(f), len(theta)))               

    # cálculo da velocidade de grupo para diferentes thetas
    for i in range(len(theta)):
        cgf[:, i] =( 1 - (c**2*(dt/dh)**2 * ((81/64)*(np.cos((2*pi*f)/c*np.sin(theta[i])*dh) + \
        np.cos((2*pi*f)/c*np.cos(theta[i])*dh) - 2) - (9/96)*(np.cos(2*(2*pi*f)/c*np.sin(theta[i])*dh)-np.cos((2*pi*f)/c*np.sin(theta[i])*dh) + \
        np.cos(2*(2*pi*f)/c*np.cos(theta[i])*dh) - np.cos((2*pi*f)/c*np.cos(theta[i])*dh))+ (1/576)*(np.cos(3*(2*pi*f)/c*np.sin(theta[i])*dh) + \
        np.cos(3*(2*pi*f)/c*np.cos(theta[i])*dh) - 2)) + 1 )**2)**(-1/2) * \
        (c**2*(dt/dh) * ((81/64)*(np.sin((2*pi*f)/c*np.sin(theta[i])*dh)*np.sin(theta[i]) + \
        np.sin((2*pi*f)/c*np.cos(theta[i])*dh) *np.cos(theta[i])) - (9/96)*(np.sin(2*(2*pi*f)/c*np.sin(theta[i])*dh)*2 *np.sin(theta[i])- \
        np.sin((2*pi*f)/c*np.sin(theta[i])*dh) *np.sin(theta[i])+ np.sin(2*(2*pi*f)/c*np.cos(theta[i])*dh)*2 *np.cos(theta[i]) - \
        np.sin((2*pi*f)/c*np.cos(theta[i])*dh) *np.cos(theta[i]))+ (1/576)*(np.sin(3*(2*pi*f)/c*np.sin(theta[i])*dh)*3 *np.sin(theta[i]) + \
        np.sin(3*(2*pi*f)/c*np.cos(theta[i])*dh)*3 *np.cos(theta[i]) )) )

    return cgf

def calc_cp(theta : np.ndarray, f : np.ndarray, dt : float, dh : float, c : float) -> np.ndarray:
    
    cpf = np.zeros((len(f), len(theta)))


    for i in range(len(theta)):
        
        cpf[:, i] = (c / (2*np.pi*f*dt)) * np.arccos(c**2*(dt/dh)**2 * ((81/64)*(np.cos((2*np.pi*f)/c*np.sin(theta[i])*dh) + \
        np.cos((2*np.pi*f)/c*np.cos(theta[i])*dh) - 2) - (9/96)*(np.cos(2*(2*np.pi*f)/c*np.sin(theta[i])*dh)- \
        np.cos((2*np.pi*f)/c*np.sin(theta[i])*dh) + np.cos(2*(2*np.pi*f)/c*np.cos(theta[i])*dh) - np.cos((2*np.pi*f)/c*np.cos(theta[i])*dh)) \
        + (1/576)*(np.cos(3*(2*np.pi*f)/c*np.sin(theta[i])*dh) + np.cos(3*(2*np.pi*f)/c*np.cos(theta[i])*dh) - 2)) + 1 )

    return cpf



def plot_cg(f : np.ndarray, cgf : np.ndarray, cre : np.ndarray, labels : list[str], 
            SaveName : str):

    plt.figure(figsize=(13, 6))

    for i, l in zip(range(len(cgf[:, 0])), labels):
        plt.plot(f, cgf[:, i], lw=2, label=l)

    plt.plot(f, cre[:, 0], ls='--', lw=3, color='black', label='real')

    plt.grid()
    plt.xlabel('k [1/m]', fontsize=14, fontweight='bold')
    plt.ylabel('C_g(k, θ) [m/s]', fontsize=14, fontweight='bold')
    plt.legend(fontsize=14)
    plt.xlim((5, 55))

    plt.tick_params(axis='both', labelsize=16)

    ax = plt.gca()

    ax.ticklabel_format(style='plain', axis='y')
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.get_major_formatter().set_useOffset(False)

    if not SaveName.endswith('.png'):
        SaveName += '.png'

    plt.savefig(SaveName, dpi=300, bbox_inches='tight')


    plt.show()



def plot_comp(x : np.ndarray, y1 : np.ndarray, y2 : np.ndarray,
              ox : int, nx : int, dx : float, labels : list[str],
              xylabel : list[str], SaveName : str):
    '''
        Plot da comparação entre dois vetores.

        Entradas:
            x  -> Eixo das abicissas
            y1 -> Primeiro dado
            y2 -> Segundo dado
            ox -> Origem do eixo
            nx -> Número de amostras do eixo
            dx -> Espaçamento das amostras do eixo
            labels -> lista com legenda dos dois dados
            xylabels -> lista com label dos eixos. (Primeiro elemento: eixo x, segundo elemento: eixo y)
        Saida:
            Gráfico de comparação entre os vetores y1 e y2
    '''
    
    plt.figure(figsize=(13, 5))

    plt.plot(x, y1, c='red', lw=2, ls='--', label=labels[0])
    plt.plot(x, y2, c='blue', lw=2, label=labels[1])

    plt.xlabel(xylabel[0], fontsize=13, fontweight='bold')
    plt.ylabel(xylabel[1], fontsize=13, fontweight='bold')

    plt.grid()

    # plt.xlim((ox, (nx * dx) - 0.75))
    plt.xlim(4, 6)

    plt.legend(fontsize=14)

    if not SaveName.endswith('.png'):
        SaveName += '.png'

    plt.savefig(SaveName, dpi=300, bbox_inches='tight')

    plt.show()

def plot(x : np.ndarray, y1 : np.ndarray,
              ox : int, nx : int, dx : float, labels : list[str],
              xylabel : list[str], SaveName : str, legend : bool = True):

    plt.figure(figsize=(13, 5))

    plt.plot(x, y1, c='black', lw=2, label=labels[0])

    plt.xlabel(xylabel[0], fontsize=13, fontweight='bold')
    plt.ylabel(xylabel[1], fontsize=13, fontweight='bold')

    plt.grid()

    # plt.xlim((ox, (nx * dx)- 0.75)) 
    plt.xlim(4, 6)

    if legend is True:
        plt.legend(fontsize=14)

    if not SaveName.endswith('.png'):
        SaveName += '.png'

    plt.savefig(SaveName, dpi=300, bbox_inches='tight')


    plt.show()


def comp_cg(f : np.ndarray, f_est : np.ndarray, cg_teo : np.ndarray, cg_est : np.ndarray,
            label1 : list[str], label2 : list[str], SaveName : str):

    n_teo = len(cg_teo[0, :])
    n_est = len(cg_est[0, :])

    plt.figure(figsize=(15,8))

    for i, l1 in zip(range(n_teo), label1):
        plt.plot(f_est, cg_est[:, i], lw=2, ls='--', label=l1)

    for j, l2 in zip(range(n_est), label2):
        plt.plot(f, cg_teo[:, j], lw=2, label=l2)

    plt.grid()
    plt.xlabel('f [Hz]', fontsize=14, fontweight='bold')
    plt.ylabel('C_g [m/s]', fontsize=14, fontweight='bold')
    plt.legend(fontsize=14)
    plt.xlim((5,55))

    plt.tick_params(axis='both', labelsize=16)

    ax = plt.gca()

    ax.ticklabel_format(style='plain', axis='y')
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.get_major_formatter().set_useOffset(False)

    if not SaveName.endswith('.png'):
        SaveName += '.png'

    plt.savefig(SaveName, dpi=300, bbox_inches='tight')

    plt.show()


