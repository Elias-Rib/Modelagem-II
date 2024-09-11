import numpy as np
import matplotlib.pyplot as plt 

#%%--------#Parâmetros Malha e Tempo#--------#
#Comprimento total do domínio
h = 0.2               #m
#Distância entre nós da malha
delta_y = h/5      #m
#Criando vetor espacial
x = np.arange(0, h+delta_y, delta_y)  
#Número total nós da malha
n = len(x)

vp = 0.5
dp = -0.3

rho = 1000    #kg/m3
mi = 0.001    #Pa.s
 
#Tempo total de simulação 
tf = 15000      #s     
#Passo de tempo
delta_t = 1     #s 
#Criando vetor temporal
t = np.arange(0, tf+delta_t, delta_t) 
#Número total de passos de tempo
p = len(t)
    
#%%--------#Preparar animação#--------#    
# Definição do gráfico
def grafico(j):
    # Gráfico colormap
    Z = np.column_stack((Vx[:, j], Vx[:, j])).T
    y = np.arange(2)
    fig, axis = plt.subplots(2, 1, figsize=[7, 5])
    plt.tight_layout()
    
    ax = axis[0]
    tg = j * delta_t
    ax.set_title('tempo [s] = %1.2f' % tg)
    cmap = plt.get_cmap('bone_r')
    ax.get_yaxis().set_visible(False)
    im = ax.pcolormesh(x, y, Z, cmap=cmap, shading='gouraud', vmin=0, vmax=max(CC1, CC2))
    clb = fig.colorbar(im)
    clb.ax.set_title('Vx (m/s)')
    
    # Gráfico dispersão
    ax = axis[1]
    ax.plot(x, Vx[:, j], 'k*')
    ax.set_ylim([0, max(CC1, CC2)])
    ax.set_xlabel('Distância, m')
    ax.set_ylabel('Vx, m/s')
    plt.pause(0.1)  # Pausar para exibição da animação
    plt.show()

#%%--------#Dados do Problema#--------#
# Condição contorno 1, em x=0 (Placa Topo)
CC1 = 1  # mol/m3 
# Condição contorno 2, em x=L (Placa Base)
CC2 = 0  # mol/m3 
# Condição inicial no tempo (todas as posições internas)
CCI = 0  # mol/m3 

#%%--------#Solver MDF#--------#
# Criar matriz resultados 
Vx = np.zeros((n, p))
# Inicializar domínio no tempo zero
Vx[:, 0] = CCI
# Condição de contorno 1
Vx[0, :] = CC1
# Condição de contorno 2
Vx[-1, :] = CC2

# Resolvendo EDP explicitamente no tempo por MDF
for j in range(1, p):  # no tempo, exceto tempo inicial
    for i in range(1, n - 1):  # no espaço, exceto contornos
        Vx[i, j] = (-dp)*delta_t/rho + mi*((Vx[i+1, j-1] - 2*Vx[i, j-1] + Vx[i-1, j-1])/delta_y**2)*delta_t/rho + Vx[i, j-1]
    
    # Chamar a função do gráfico para animar
    if j % 100 == 0:  # Para atualizar o gráfico a cada 100 passos de tempo
        grafico(j)

plt.figure(1)
plt.plot(Vx[:,-1],x)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('velocidade, m/s')
plt.ylabel('altura, m')
plt.title('Escoamento Couette')