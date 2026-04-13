import numpy as np
import util
import os
import matplotlib.pyplot as plt

SaveDir = '../figs_comp'
SavePath = os.path.join(os.getcwd(), SaveDir)

if not os.path.exists(SavePath):
    os.makedirs(SavePath)

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

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(13, 5))

extent = [x.min(), x.max(), z.max(), z.min()]

im = plt.imshow(vel, cmap='viridis', extent=extent, vmin=900, vmax=1000, aspect=0.5)

# Color bar
cbar = plt.colorbar(im, shrink=1, pad=0.05)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Velocity (m/s)', fontsize=13, fontweight='bold')

# Fonte
plt.scatter(fonte[0], fonte[1], color='red', marker='*', s=150, label='Source')

# Receptores
plt.scatter(rec_1[0, 0], rec_1[0, 1], color='blue', marker='^', s=100, label='θ = 0')    # primeira linha
plt.scatter(rec_2[0, 0], rec_2[0, 1], color='blue', marker='^', s=100,)                     # segunda linha

plt.scatter(rec_1[1, 0], rec_1[1, 1], color='black', marker='^', s=100, label='θ = π/12')    # primeira linha
plt.scatter(rec_2[1, 0], rec_2[1, 1], color='black', marker='^', s=100,)                     # segunda linha

plt.scatter(rec_1[2, 0], rec_1[2, 1], color='green', marker='^', s=100, label='θ = π/6')    # primeira linha
plt.scatter(rec_2[2, 0], rec_2[2, 1], color='green', marker='^', s=100,)                     # segunda linha

plt.scatter(rec_1[3, 0], rec_1[3, 1], color='purple', marker='^', s=100, label='θ = π/4')    # primeira linha
plt.scatter(rec_2[3, 0], rec_2[3, 1], color='purple', marker='^', s=100,)                     # segunda linha


# Labels
plt.xlabel('Distance (m)', fontsize=13, fontweight='bold', labelpad=10)
plt.ylabel('Depth (m)', fontsize=13, fontweight='bold', labelpad=20)

# Eixos
plt.tick_params(axis='both', labelsize=14)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')

plt.legend(fontsize=12, loc='upper left')

# Salvar figura
plt.savefig(os.path.join(SavePath, 'modelo_cg_cp.png'), dpi=300, bbox_inches='tight')
plt.show()

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

frq = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]   # frequências de medida
dr  = 300                                           # distância entre os receptores


cre = np.ones((len(f), 1))                          # vetor para velocidade real para f
cre *= 1000                                         # vetor com a velocidade real do meio


cgf = util.calc_cg(theta, f, dt, dh, c)

util.plot_cg(
    f,
    cgf, cre,
    ['θ = 0', 'θ = π/12', 'θ = π/6', 'θ = π/4'],
    os.path.join(SavePath, 'group_vel')
)

# salva os dados
for i in range(0, len(theta)):
    np.savetxt(
        f"dados/cgf_{i}.txt",
        np.column_stack((f, cgf[:, i])),
        fmt="%.3f"
    )

# ==========================
#      Medição de Cg
# ==========================

data1  = []
data2  = []
cg_est = []

DirFrq = os.path.join(os.getcwd(), '../dados_frq')


# TODO: ATENÇÃO A ORDEM DOS LOOPS VERIFICAR SE MUDA O RESULTADO SE INVERTER.

for f in frq:
    for i in range(0, len(theta)):
        data_mod1 = util.loadrsf(os.path.join(DirFrq, f'dataP1-{(i+1)}_fos{f}'))

        data1.append(data_mod1)

    for i in range(4, 2*len(theta)):
        data_mod2 = util.loadrsf(os.path.join(DirFrq, f'dataP2-{(i+1)}_fos{f}'))

        data2.append(data_mod2)

data1 = np.array(data1).reshape(len(frq), len(theta), nt)
data2 = np.array(data2).reshape(len(frq), len(theta), nt)



for d1, d2 in zip(data1, data2):
    for i in range(0, len(theta)):

        # Correlação cruzada entre os dados
        corr = np.correlate(d2[i, :], d1[i, :], mode='full')

        # lags para criar eixo do tempo
        lags = np.arange(-len(d1[i, :])+1, len(d1[i, :]))

        # eixo do tempo
        t_corr = lags * dt


        # util.plot(
        #     t_corr, corr, 0, len(t_corr), (t_corr[1] - t_corr[0]), ['Cross Correlation'],
        #     ['time_lag (s)', 'Cross-correlation'], os.path.join(SavePath, f'corr_{i}')
        # )

        # indicie do pico máximo da correlação
        idx = np.argmax(corr)

        max_t  = lags[idx] * dt
        max_t_ = lags[idx-1] * dt
        max_t1 = lags[idx+1] * dt


        # plt.figure(figsize=(10,5))
        # plt.plot(t_corr, corr, c='black', lw=2)

        # plt.scatter(max_t, corr[idx], c='red', marker="o")
        # plt.scatter(max_t_, corr[idx-1], c='red', marker="o")
        # plt.scatter(max_t1, corr[idx+1], c='red', marker="o")        
        
        # plt.grid()
        # plt.xlim((0.29, 0.31))
        # plt.show()

        if 0 < idx < len(corr)-1:
            y1 = corr[idx-1]
            y2 = corr[idx]
            y3 = corr[idx+1]

            delta = 0.5 * ((y1 - y3) / (y1 - 2*y2 + y3))
        else:
            delta = 0.0

        
        lag_r = lags[idx] + delta
        delta_t = lag_r * dt

        # velocidade de grupo medida
        cg_med = dr / delta_t

        # adiciona valor ao vetor
        cg_est.append(cg_med)


cg_est  = np.array(cg_est).reshape(len(frq), len(theta), 1)
frq     = np.array(frq).reshape(len(frq), 1)

for i in range(0, len(theta)):
    # salva valores medidos
    np.savetxt(f'dados/dados_cg_{i}.txt', 
        np.column_stack((frq, cg_est[:, i])),
        fmt='%.3f')
    

# Seção de plot e salvamento

# Carregar dados
f_teo, cg_teo = util.load_data_cg(theta, 'cgf')
f_est, cg_est = util.load_data_cg(theta, 'dados_cg')

util.comp_cg(
    f_teo, f_est,
    cg_teo, cg_est,
    ['Measured (θ = 0)', 'Measured (θ = π/12)', 'Measured (θ = π/6)', 'Measured (θ = π/4)'],
    ['Theoretical (θ = 0)', 'Theoretical (θ = π/12)', 'Theoretical (θ = π/6)', 'Theoretical (θ = π/4)'],
    os.path.join(SavePath, 'cg_teo_med')
)