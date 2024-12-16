import numpy as np
import matplotlib.pyplot as plt

# Propriedades das camadas
camadas = [
    {"h": 0.0001, "k": 0.235, "rho": 1200, "c": 3500, "omega_b": 0.0, "Q_m": 0.0, "q": 2.0},  # Epiderme
    {"h": 0.0007, "k": 0.445, "rho": 1000, "c": 3000, "omega_b": 0.0002, "Q_m": 368.1, "q": 2.0},  # Derme
    {"h": 0.0008, "k": 0.445, "rho": 900,  "c": 2500, "omega_b": 0.0013, "Q_m": 368.1, "q": 2.0},  # Gordura
    {"h": 0.002, "k": 0.185, "rho": 1050, "c": 3600, "omega_b": 0.0001, "Q_m": 368.3, "q": 2.0},  # Músculo superficial
    {"h": 0.005, "k": 0.558, "rho": 1030, "c": 3852, "omega_b": 0.0063, "Q_m": 3700, "q": 1.1},
    {"h": 0.008, "k": 0.51, "rho": 1080, "c": 3700, "omega_b": 0.0027, "Q_m": 684.2, "q": 2.0},  # Músculo profundo
]

# Parâmetros gerais
T_a = 37  # Temperatura arterial [°C]
T_ext = 25  # Temperatura externa [°C]
rho_b = 1060
c_pb = 4200 
h_ext = 100  

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
q =  np.zeros_like(T)

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
    q[indices] = camada["q"]

# Funções para cálculo de Q e w
def calc_Q(T, Q_bas, q):
    """Cálculo da taxa metabólica ajustada"""
    return Q_bas * q ** ((T - T_a) / 10)

def calc_w(T, omega, q):
    """Cálculo da perfusão sanguínea ajustada"""
    return omega * q ** ((T - T_a) / 10)
w_m = calc_w(T_a, omega_b, q)
Q_m = calc_Q(T_a, Q_m_values, q)

# Listar tempos para plotar
tempos = [0, int(nt/4), int(nt/2), int(3*nt/4), nt-1]
temperaturas = []

# Simulação temporal
for n in range(nt):
    T_new = T.copy()  # Criação de uma nova variável para armazenar as temperaturas futuras
        # Condição de contorno à esquerda (Eq. 5: fluxo convectivo)
    T_new[0] = T_ext

    # Condição de contorno à direita (temperatura fixada em T_a)
    T_new[-1] = T_a
    
    for i in range(1, nx_total - 1):  # Para cada ponto da malha, exceto os limites
        # Termo de condução térmica (difusividade)
        conduction = (k[i] * (T[i + 1] - 2 * T[i] + T[i - 1])) / dx**2
        # Termo de convecção (efetivo da perfusão sanguínea)
        convection = w_m[i] * rho_b * c_pb * (T_a - T[i])
        # Fonte de calor (Q_m)
        heat_source = Q_m[i]
        
        # Atualização da temperatura no ponto i usando a equação discretizada
        T_new[i] = T[i] + dt * (conduction + convection + heat_source) / (rho[i] * c[i])

    # Condições nas interfaces (Eq. 6)
    for j in range(1, len(camadas)):
        interface_idx = int(sum(camada["h"] for camada in camadas[:j]) / comprimento_total * nx_total)
        T_new[interface_idx] = (k[interface_idx] * T[interface_idx + 1] + k[interface_idx - 1] * T[interface_idx - 1]) / (k[interface_idx] + k[interface_idx - 1])

    # Atualiza a temperatura com os novos valores
    T = T_new.copy()

    # Armazenar temperaturas para os tempos especificados
    if n in tempos:
        temperaturas.append(T.copy())

# Visualização
plt.figure(figsize=(6, 6))
for i, temp in enumerate(temperaturas):
    plt.plot(x, temp, label=f'Tempo = {tempos[i]*dt:.2f} s')

plt.xlabel('Distância [m]')
plt.ylabel('Temperatura [°C]')
plt.title('Distribuição de Temperatura ao Longo do Tempo')
plt.legend()
plt.grid(True)
plt.show()
