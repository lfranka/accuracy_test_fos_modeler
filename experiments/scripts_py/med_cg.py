import numpy as np
import util
import os
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.ndimage import gaussian_filter1d



SaveDir = '../figs'
SavePath = os.path.join(os.getcwd(), SaveDir)

if not os.path.exists(SavePath):
    os.makedirs(SavePath)


Datadir     = 'dados'
DataPath    = os.path.join(os.getcwd(), Datadir)

if not os.path.exists(DataPath):
    os.makedirs(DataPath)


pi = np.pi

# ============================================
#      Plotando o modelo de velocidade
# ============================================
# Parâmetros do modelo
nx, nz = 701, 701
dx, dz = 2.0, 2.0

theta = [0, np.pi/12, np.pi/6, np.pi/4]


# Velocidade do modelo
vel = np.ones((nz, nx)) * 1000  # m/s

# Criando os eixos x e z
x = np.linspace(0, (nx-1)*dx, nx)
z = np.linspace(0, (nz-1)*dz, nz)

# Posição da fonte e receptor  
fonte    = [300, 500]

rec_1 = np.zeros((len(theta), 2))   # primeira linha de receptores
rec_2 = np.zeros((len(theta), 2))   # segunda linha de receptores

r   = 700.2     # m
dr  = 300       # m

for i in range(0, len(theta)):
    
    rx = fonte[0] + (r * np.cos(theta[i]))
    rz = fonte[1] + (r * np.sin(theta[i]))

    rx_ = fonte[0] + ((r + dr) * np.cos(theta[i]))
    rz_ = fonte[1] + ((r + dr) * np.sin(theta[i]))

    rec_1[i, 0] = rx
    rec_1[i, 1] = rz

    rec_2[i, 0] = rx_
    rec_2[i, 1] = rz_

# # ----------------------------
# # Plot
# # ----------------------------
# plt.figure(figsize=(13, 5))

# extent = [x.min(), x.max(), z.max(), z.min()]

# im = plt.imshow(vel, cmap='viridis', extent=extent, vmin=900, vmax=1000, aspect=0.5)

# # Color bar
# cbar = plt.colorbar(im, shrink=1, pad=0.05)
# cbar.ax.tick_params(labelsize=14)
# cbar.set_label('Velocity (m/s)', fontsize=13, fontweight='bold')

# # Fonte
# plt.scatter(fonte[0], fonte[1], color='red', marker='*', s=150, label='Source')

# # Receptores
# plt.scatter(rec_1[0, 0], rec_1[0, 1], color='blue', marker='^', s=100, label='θ = 0')    # primeira linha
# plt.scatter(rec_2[0, 0], rec_2[0, 1], color='blue', marker='^', s=100,)                     # segunda linha

# plt.scatter(rec_1[1, 0], rec_1[1, 1], color='black', marker='^', s=100, label='θ = π/12')    # primeira linha
# plt.scatter(rec_2[1, 0], rec_2[1, 1], color='black', marker='^', s=100,)                     # segunda linha

# plt.scatter(rec_1[2, 0], rec_1[2, 1], color='green', marker='^', s=100, label='θ = π/6')    # primeira linha
# plt.scatter(rec_2[2, 0], rec_2[2, 1], color='green', marker='^', s=100,)                     # segunda linha

# plt.scatter(rec_1[3, 0], rec_1[3, 1], color='purple', marker='^', s=100, label='θ = π/4')    # primeira linha
# plt.scatter(rec_2[3, 0], rec_2[3, 1], color='purple', marker='^', s=100,)                     # segunda linha


# # Labels
# plt.xlabel('Distance (m)', fontsize=13, fontweight='bold', labelpad=10)
# plt.ylabel('Depth (m)', fontsize=13, fontweight='bold', labelpad=20)

# # Eixos
# plt.tick_params(axis='both', labelsize=14)
# plt.gca().xaxis.tick_top()
# plt.gca().xaxis.set_label_position('top')

# plt.legend(fontsize=12, loc='upper left')

# # Salvar figura
# plt.savefig(os.path.join(SavePath, 'modelo_cg_cp.png'), dpi=300, bbox_inches='tight')
# plt.show()

# ==============================================
#   Calculando e plotando a velocidade de grupo
# ==============================================
f = np.arange(5, 55 + 0.5, 0.5)                     # frequências
dh = 2                                              # discretização espacial
dt = 0.0004                                         # discretização temporal
c = 1000                                            # velocidade real do meio
r = 700.2                                           # distância fonte receptor
theta = [0, np.pi/12, np.pi/6, np.pi/4]             # angulos entre fonte receptor
nt  = 5000
rc = r / c
rho = 1200.0

frq = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]   # frequências de medida
dr  = 300                                           # distância entre os receptores


cre = np.ones((len(f), 1))                          # vetor para velocidade real para f
cre *= 1000                                         # vetor com a velocidade real do meio


cgf = util.calc_cg(theta, f, dt, dh, c)

# util.plot_cg(
#     f,
#     cgf, cre,
#     ['θ = 0', 'θ = π/12', 'θ = π/6', 'θ = π/4'],
#     os.path.join(SavePath, 'group_vel')
# )

# salva os dados
for i in range(0, len(theta)):
    np.savetxt(
        os.path.join(DataPath,f'cgf_{i}.txt'),
        np.column_stack((f, cgf[:, i])),
        fmt="%.3f"
    )


# Medição da velocidade de grupo

t = np.arange(nt) * dt

data_mod1 = []
data_mod2 = []
data_ana  = []

cg_est   = []
del_t_list = []
qualidade_list = []

for fe in frq:
    for th in range(1, len(theta) + 1):
        data1 = util.loadrsf(f'../dados_frq/dataP1-{th}_fos{fe}.rsf')

        data_mod1.append(data1)
    
    for th in range(5, 2*len(theta)+1):
        data2 = util.loadrsf(f'../dados_frq/dataP2-{th}_fos{fe}.rsf')

        data_mod2.append(data2)


data_mod1 = np.array(data_mod1).reshape(len(frq), len(theta), len(t)) 
data_mod2 = np.array(data_mod2).reshape(len(frq), len(theta), len(t)) 


for fe in frq:
    wav = util.loadrsf(f'../dados_frq/wav_{fe}.rsf')

    # Criando a função de Green da equação de segunda ordem
    G = np.zeros_like(t)

    mask = t > rc
    G[mask] = ((-1)/(2*pi)) * (1 / np.sqrt(t[mask]**2 - rc**2))

    # Criar dado semi-analítico
    Sap = util.create_sap(G, wav, t, rho, dx, dz)

    data_ana.append(Sap)


data_ana = np.repeat(np.array(data_ana)[:, np.newaxis, :], len(theta), axis=1)


for i in range(0, len(frq)):
    for j in range(0, len(theta)):

        S1 = data_mod1[i, j, :]
        S2 = data_mod2[i, j, :]

        fft_s1 = np.fft.fft(S1)
        fft_s2 = np.fft.fft(S2)

        cross = fft_s2 * np.conj(fft_s1)        # Correlação cruzada na frequência

        n_f = np.fft.fftfreq(len(S1), dt)   # Eixo da frequências


        f_pos = n_f > 0     # Somente frequências positivas

        # Pegando somente frequências positivas e em torno da frequência central
        f_min = frq[i] * 0.8
        f_max = frq[i] * 1.2

        f_band = (n_f >= f_min) & (n_f <= f_max)

        mask = f_pos & f_band

        # supondo onda plana isso deveria ser -2*pi*f*δt + b
        phi_cross = np.unwrap(np.angle(cross[mask]))
        freqs = n_f[mask] 

        # Aproximação por mínimos quadrados
        # Aproximação de primeiro grau (estima a e b)  nesse caso a = -2*pi*δt
        a = np.polyfit(freqs, phi_cross, 1)

        # logo δt = - a / 2 * pi
        delay = -a[0] / (2 * pi)

        if  i >= 2: 
            del_t = delay
        else:   # Usa o método temporal nas primeira frequências. Continua sendo mais preciso !!!!
            corr = np.correlate(S2, S1, mode='full')
            lags = np.arange(-len(S2)+1, len(S2))
            idt = np.argmax(corr)
            
            if 0 < idt < len(corr) - 1:
                y1, y2, y3 = corr[idt-1], corr[idt], corr[idt+1]
                delta = 0.5 * ((y1 - y3) / (y1 - 2*y2 + y3))
                lag_r = lags[idt] + delta
                del_t = lag_r * dt
            else:
                del_t = lags[idt] * dt
            
        
        # Calcula velocidade de grupo
        cg_mod = dr / del_t
        cg_est.append(cg_mod)


cg_est = np.array(cg_est).reshape(len(frq), len(theta))

util.comp_cg(
    f, frq,
    cgf, cg_est,
    ['Measured (θ = 0)', 'Measured (θ = π/12)', 'Measured (θ = π/6)', 'Measured (θ = π/4)'],
    ['Theoretical (θ = 0)', 'Theoretical (θ = π/12)', 'Theoretical (θ = π/6)', 'Theoretical (θ = π/4)'],
    os.path.join(SavePath, 'cg_teo_med')
)