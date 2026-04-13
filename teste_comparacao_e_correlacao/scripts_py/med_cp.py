import numpy as np
import os
import glob
import util

SaveDir = '../figs_comp'
SavePath = os.path.join(os.getcwd(), SaveDir)

if not os.path.exists(SavePath):
    os.makedirs(SavePath)


dh = 2
dt = 0.0004
c = 1000
r = 700.2
theta = [0, np.pi/12, np.pi/6, np.pi/4]

f = np.arange(5, 55 + 0.5, 0.5)

cre = np.ones((len(f), 1))

cre *= 1000


cpf = util.calc_cp(theta, f, dt, dh, c)


util.plot_cg(
    f, cpf, cre,
    ['θ = 0', 'θ = π/12', 'θ = π/6', 'θ = π/4'],
    os.path.join(SavePath, 'pahse_vel')
)

for i in range(0,len((theta))):
    np.savetxt(
        f"dados/cpf_{i}.txt",
        np.column_stack((f, cpf[:, i])),
        fmt="%.3f"
    )



# =========================
# parâmetros
# =========================

frq = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
dr = 300   
nt = 5000

# =========================
# leitura dos dados
# =========================

data1 = []
data2 = []

DirFrq = os.path.join(os.getcwd(), '../dados_frq')

for f in frq:
    for i in range(0, len(theta)):
        data_mod1 = util.loadrsf(os.path.join(DirFrq, f'dataP1-{(i+1)}_fos{f}'))

        data1.append(data_mod1)

    for i in range(4, 2*len(theta)):
        data_mod2 = util.loadrsf(os.path.join(DirFrq, f'dataP2-{(i+1)}_fos{f}'))

        data2.append(data_mod2)


data1 = np.array(data1).reshape(len(frq), len(theta), nt)
data2 = np.array(data2).reshape(len(frq), len(theta), nt)

# Criar eixo das frequências e calcular omega
n_f = np.fft.fftfreq(nt, d=dt)
omega = 2 * np.pi * n_f

# =========================
# FFT
# =========================
for d1, d2 in zip(data1, data2):
    for i in range(0, len(theta)):

        fft1 = np.fft.fft(d1)
        fft2 = np.fft.fft(d2)


        # =========================
        # cross-spectrum correto
        # =========================
        cross = fft2 * np.conj(fft1)

        # =========================
        # fase
        # =========================

        phase_cross = np.angle(cross)
        phase_unw = np.unwrap(phase_cross)

        # =========================================================
        # medindo a velocidade de fase
        # =========================================================
        for f in frq:
            if f == frq[0]:
                f_min = f
            else:
                f_min = f + 0.5
            
            f_max = f + 5

            mask = (n_f >= f_min) & (n_f <= f_max)


            freq_band   = n_f[mask]
            omega_band  = omega[mask]
            dphi_band   = phase_unw[i, mask]

            k_num = dphi_band / dr

            cp_med = np.abs(omega_band / k_num)

            
            np.savetxt(
                f"dados/cp-{i+1}_{f}hz.txt",
                np.column_stack((freq_band, cp_med)),
                fmt="%.3f"
            )


f_list = []
cp_list = []
for i in range(0, len(theta)):
    for f in frq:
        file = np.loadtxt(f'dados/cp-{i+1}_{f}hz.txt', dtype='f')

        f_list.append(file[:,0])
        cp_list.append(file[:,1])


f_total = np.concatenate(f_list).reshape(4, 101, 1)
cp_total = np.concatenate(cp_list).reshape(4, 101, 1)

# Remover arquivos
files = glob.glob('dados/cp-*hz.txt')
for f in files:
    os.remove(f)

for i in range(0, len(theta)):
    np.savetxt(
        f"dados/dados_cp_{i}.txt",
        np.column_stack((f_total[i, :], cp_total[i, :])),
        fmt="%.3f"
    )


# Plot e salvamento dos dados
f_teo, cg_teo = util.load_data_cg(theta, 'cpf')
f_est, cg_est = util.load_data_cg(theta, 'dados_cp')

util.comp_cg(
    f_teo, f_est,
    cg_teo, cg_est,
    ['Measured (θ = 0)', 'Measured (θ = π/12)', 'Measured (θ = π/6)', 'Measured (θ = π/4)'],
    ['Theoretical (θ = 0)', 'Theoretical (θ = π/12)', 'Theoretical (θ = π/6)', 'Theoretical (θ = π/4)'],
    os.path.join(SavePath, 'cp_teo_med')
)