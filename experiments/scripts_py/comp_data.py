import numpy as np
import os
import util
import matplotlib.pyplot as plt

SaveDir = '../figs'
SavePath = os.path.join(os.getcwd(), SaveDir)

if not os.path.exists(SavePath):
    os.makedirs(SavePath)

# ============================================
#      Plotando o modelo de velocidade
# ============================================
# Parâmetros do modelo
nx, nz = 701, 701
dx, dz = 2.0, 2.0

# Velocidade do modelo
vel = np.ones((nz, nx)) * 1000  # m/s

# Criando os eixos x e z
x = np.linspace(0, (nx-1)*dx, nx)
z = np.linspace(0, (nz-1)*dz, nz)

# Posição da fonte e receptor  
fonte = (300, 500)          # m
receptor = (1000.2, 500)    # m

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

# Receptor
plt.scatter(receptor[0], receptor[1], color='blue', marker='^', s=100, label='Receiver')

plt.xlabel('Distance (m)', fontsize=13, fontweight='bold', labelpad=10)
plt.ylabel('Depth (m)', fontsize=13, fontweight='bold', labelpad=20)


plt.tick_params(axis='both', labelsize=14)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')

plt.legend(fontsize=12, loc='upper left')

plt.savefig(os.path.join(SavePath, 'modelo_comp.png'), dpi=300, bbox_inches='tight')
plt.show()


# ============================================
#                Experimentos
# ============================================
# Parâmetros
r    = 700.2
r_rc = 700.4
r_dh = 700.25
c    = 1000.0
rho  = 1200.0
pi   = np.pi

frq = [5, 10, 20, 45]

rc      = r / c
rc_rc   = r_rc / c
rc_dh   = r_dh / c

dx = dz = 2

nt = 5000
ot = 0
dt = 0.0004

# Criar eixo do tempo
t = np.arange(nt) * dt

# Criando a função de Green da equação de segunda ordem
G = np.zeros_like(t)

mask = t > rc
G[mask] = ((-1)/(2*pi)) * (1 / np.sqrt(t[mask]**2 - rc**2))

# ============================================
#           Comparação fos e fo
# ============================================

# Diretórios com os arquivos .rsf
DirFo   = os.path.join(os.getcwd(), '../dados_fo')

# Retirar os dados
data_fo     = util.loadrsf(os.path.join(DirFo, 'dataP_fo'))
data_fos    = util.loadrsf(os.path.join(DirFo, 'dataP_fos'))
wav_fo      = util.loadrsf(os.path.join(DirFo, 'wav_fo'))

# Criar sinal semi-analítico
Sap_fo_fos = util.create_sap(G, wav_fo, t, rho, dx, dz)

# Plot das comparações
util.plot_comp(
    t, data_fo, Sap_fo_fos,
    ot, nt, dt, ['Modeled (fo)', 'Semi-analytical'],
    ['t (s)', 'Amplitude'], os.path.join(SavePath, 'comp_fo')
)

util.plot_comp(
    t, data_fos, Sap_fo_fos,
    ot, nt, dt, ['Modeled (fos)', 'Semi-analytical'],
    ['t (s)', 'Amplitude'], os.path.join(SavePath, 'comp_fos')
)



# ============================================
#    Comparação modelado X semi-analítico
# ============================================
# Erro relativo
erro = util.diff_data(Sap_fo_fos, data_fos)

util.plot(
    t, erro, ot, nt, dt,
    ['Reltativo Error'], ['t (s)', 'Error (%)'], 
    os.path.join(SavePath, 'diff_fos')
)

# ============================================
#       Localização da Função de Green
# ============================================
G_rc = np.zeros_like(t)

mask = t > rc_rc
G_rc[mask] = ((-1)/(2*pi)) * (1 / np.sqrt(t[mask]**2 - rc_rc**2))

# Diretório com os arquivos .rsf
DirRc = os.path.join(os.getcwd(), '../dados_rc')

# Retirar os dados
data_rc = util.loadrsf(os.path.join(DirRc, 'dadoP_rc'))
wav_rc  = util.loadrsf(os.path.join(DirRc, 'wav_rc'))

# Criar sinal semi-analítico
Sap_rc = util.create_sap(G_rc, wav_rc, t, rho, dx, dz)

# plot dos resultados
util.plot_comp(
    t, data_rc, Sap_rc,
    ot, nt, dt, ['Modeled (rc)', 'Semi-analitycal'],
    ['t (s)', 'Amplitude'], os.path.join(SavePath, 'comp_rc')
)

# Erro realtivo
util.plot(
    t, util.diff_data(Sap_rc, data_rc), ot, nt, dt,
    ['Relative Error'], ['t (s)', 'Error (%)'],
    os.path.join(SavePath, 'diff_rc')
)

# =============================================
#          Discretização do modelo
# =============================================
nt_dh = 20000
ot_dh = 0
dt_dh = 0.0001

dx_dh = 1
dz_dh = 1

t_dh  = np.arange(nt_dh) * dt_dh

# Criando a função de Green da equação de segunda ordem
G_dh = np.zeros_like(t_dh)

mask = t_dh > rc_dh
G_dh[mask] = ((-1)/(2*pi)) * (1 / np.sqrt(t_dh[mask]**2 - rc_dh**2))

# Diretório com arquivos .rsf
DirDh = os.path.join(os.getcwd(), '../dados_dh')

# Retirar dados
data_dh = util.loadrsf(os.path.join(DirDh, 'dataP_dh'))
wav_dh  = util.loadrsf(os.path.join(DirDh, 'wav_dh'))

# Criar dado semi-analítico
Sap_dh = util.create_sap(G_dh, wav_dh, t_dh, rho, dx_dh, dz_dh)

# Plot dos resultados
util.plot_comp(
    t_dh, data_dh, Sap_dh,
    ot_dh, nt_dh, dt_dh, ['Modeled (dh)', 'Semi-analytical'],
    ['t (s)', 'Amplitude'], os.path.join(SavePath, 'comp_dh')
)

util.plot(
    t_dh, util.diff_data(Sap_dh, data_dh), 
    ot_dh, nt_dh, dt_dh,
    ['Relative Error'], ['t (s)', 'Error (%)'],
    os.path.join(SaveDir, 'diff_dh')
)

# =============================================
#       Variação da frequência central
# =============================================

# Diretório com arquivos .rsf
DirFrq = os.path.join(os.getcwd(), '../dados_frq')

for f in frq:
    # Retirar dados
    data_frq = util.loadrsf(os.path.join(DirFrq, f'dataP1-1_fos{f}'))
    wav_frq  = util.loadrsf(os.path.join(DirFrq, f'wav_{f}'))

    # Criar dado semi-analítico
    Sap_frq = util.create_sap(G, wav_frq, t, rho, dx, dz)

    # plot dos resultados
    util.plot_comp(
        t, data_frq, Sap_frq,
        ot, nt, dt, [f'Modeled ({f} Hz)', 'Semi-analitycal'],
        ['t (s)', 'Amplitude'], os.path.join(SavePath, f'comp_{f}hz')
    )

    # Erro realtivo
    util.plot(
        t, util.diff_data(Sap_frq, data_frq), ot, nt, dt,
        ['Relative Error'], ['t (s)', 'Error (%)'],
        os.path.join(SavePath, f'diff_{f}hz')
    )
