import numpy as np
import matplotlib.pyplot as plt

# Propriedades das camadas
camadas = [
    {"h": 0.008, "k": 0.60, "rho": 1080, "c": 3700, "omega_b": 0.0027, "Q_m": 684.2},  # Músculo profundo
    {"h": 0.002, "k": 0.50, "rho": 1050, "c": 3600, "omega_b": 0.0001, "Q_m": 368.3},  # Músculo superficial
    {"h": 0.005, "k": 0.558, "rho": 1030, "c": 3852, "omega_b": 0.0063, "Q_m": 3700},  # Tumor 
    {"h": 0.0008, "k": 0.15, "rho": 900,  "c": 2500, "omega_b": 0.0013, "Q_m": 368.1},  # Gordura
    {"h": 0.0007, "k": 0.30, "rho": 1000, "c": 3000, "omega_b": 0.0002, "Q_m": 368.1},  # Derme
    {"h": 0.0001, "k": 0.45, "rho": 1200, "c": 3500, "omega_b": 0.0, "Q_m": 0.0},  # Epiderme
]

# Parâmetros gerais
T_a = 37  # Temperatura arterial [°C]
T_ext = 25 # Temperatura externa [°C]

# Configuração da malha e do domínio
nx_total = 100  # Total de pontos na malha
comprimento_total = sum(camada["h"] for camada in camadas)
x = np.linspace(0, comprimento_total, nx_total)
dx = x[1] - x[0]
dt = 0.01  # Passo de tempo inicial (será corrigido)
nt = 10000  # Número de passos de tempo

# Inicializar matriz de temperatura
T = np.linspace(20, T_a, nx_total)

# Criar mapeamento de propriedades
k = np.zeros_like(T)
rho = np.zeros_like(T)
c = np.zeros_like(T)
omega_b = np.zeros_like(T)
Q_m_values = np.zeros_like(T)

# Atribuindo as propriedades de cada camada
x_position = 0
for camada in camadas:
    indices = (x >= x_position) & (x < x_position + camada["h"])
    k[indices] = camada["k"]
    rho[indices] = camada["rho"]
    c[indices] = camada["c"]
    omega_b[indices] = camada["omega_b"]
    Q_m_values[indices] = camada["Q_m"]
    x_position += camada["h"]

# Listar tempos para plotar
tempos = [0, int(nt/4), int(nt/2), int(3*nt/4), nt-1]
temperaturas = []

# Simulação temporal
for n in range(nt):
    T_new = T.copy()  # Criação de uma nova variável para armazenar as temperaturas futuras
    
    for i in range(1, nx_total - 1):  # Para cada ponto da malha, exceto os limites
        # Termo de condução térmica (difusividade)
        conduction = k[i] * (T[i + 1] - 2 * T[i] + T[i - 1]) / dx**2
        # Termo de convecção (efetivo da perfusão sanguínea)
        convection = omega_b[i] * (T_a - T[i])
        # Fonte de calor (Q_m)
        heat_source = Q_m_values[i]
        
        # Atualização da temperatura no ponto i usando a equação discretizada
        T_new[i] = T[i] + dt * (conduction + convection + heat_source) / (rho[i] * c[i])

    # Condição de contorno à esquerda (fluxo de calor, com temperatura externa T_ext)
    T_new[0] = T_ext
    
    # Condição de contorno à direita (temperatura fixada em T_a)
    T_new[-1] = T_a

    # Atualiza a temperatura com os novos valores
    T = T_new.copy()

    # Armazenar temperaturas para os tempos especificados
    if n in tempos:
        temperaturas.append(T.copy())

# Visualização
plt.figure(figsize=(10, 6))
for i, temp in enumerate(temperaturas):
    plt.plot(x, temp, label=f'Tempo = {tempos[i]*dt:.2f} s')

plt.xlabel('Distância [m]')
plt.ylabel('Temperatura [°C]')
plt.title('Distribuição de Temperatura ao Longo do Tempo')
plt.legend()
plt.grid(True)
plt.show()