import math
import sys
import matplotlib.pyplot as plt

"""
Cálculo de F(x)
Recebe o valor de x, float
Retorna o valor de F(x), float
"""
def F(x):
    return math.sin(x**2) + 1.1 - math.exp(-x)

"""
Cálculo do valor de F'(x)
Recebe o valor de x, float
Retorna o valor de F'(x), float
"""
def F_derivada(x):
    return 2*x*math.cos(x**2) + math.exp(-x)

"""
Aplicação do método de bisseção
Recebe valores de a, b e a tolerância eps, todos float
Retorna o valor da raiz aproximado (float), valor do erro majorado (float), número de iterações (int)
"""
def Bissecao(a, b, eps):
    fa = F(a)
    fb = F(b)   
    dfa = F_derivada(a)
    dfb = F_derivada(b)
    if fa*fb >= 0 and dfa*dfb <= 0:
        raise ValueError("F(a) e F(b) devem ter sinais opostos e as suas derivadas o mesmo sinal. Método da bisseção abortado.")   
    else:
        print(" F(a) e F(b) têm sinais opostos, F'(a) e F'(b) têm o mesmo sinal — condições verificadas para a bisseção.")
    iteracoes = 0
    while (b-a)/2 > eps:
        c = (a+b)/2
        fc = F(c)
        if fc == 0:
            return c, 0, iteracoes  
        elif fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        iteracoes += 1
    x_aproximado = (a+b)/2
    erro_majorado = (b-a)/2
    return x_aproximado, erro_majorado, iteracoes

"""
Cálculo do majorante do erro
Recebe os valores de a e b (float)
Retorna o valor do majorante do erro (float)
"""
def Majorante_erro(a, b, samples = 2000):
    espacamento = (b-a)/samples
    m_erro = float('inf')
    x = a
    for _ in range(samples + 1):
        m_erro = min(m_erro, abs(F_derivada(x)))
        x += espacamento
    if m_erro < 1e-14:
        m_erro = 1e-14
    return m_erro

"""
Aplicação do método de Newton
Recebe os valores de a, b, eps, todos float
Retorna os valor aproximado da raiz (float), o valor do erro majorado (float), número de iterações (int)
"""
def Newton(a, b, eps, x0=None):
    fa = F(a)
    fb = F(b)   
    dfa = F_derivada(a)
    dfb = F_derivada(b)
    if fa*fb >= 0 and dfa*dfb <= 0:
        raise ValueError("F(a) e F(b) devem ter sinais opostos e as suas derivadas o mesmo sinal. Método de Newton abortado.")   
    else:
        print(" F(a) e F(b) têm sinais opostos, F'(a) e F'(b) têm o mesmo sinal — condições verificadas para Newton.")
    if x0 is None:
        x0 = (a+b)/2
    m_erro = Majorante_erro(a, b)
    x = x0
    max_iteracoes = 15
    iteracoes = 0    
    while iteracoes < max_iteracoes:
        fx = F(x)
        dfx = F_derivada(x)
        if abs(dfx) < 1e-14:
            x1 = (a+b)/2
        else:
            x1 = x - (fx/dfx)
        if not (a <= x1 <= b):
            x1 = (a+b)/2
        iteracoes += 1
        erro_majorado = abs(F(x1)) / m_erro
        if erro_majorado <= eps:
            return x1, erro_majorado, iteracoes
        x = x1
    return x, abs(F(x))/m_erro, iteracoes

"""
Plota o gráfico de F(x) e encontra automaticamente o primeiro intervalo onde há mudança de sinal
Recebe xmin, xmax, espacamento, todos float
Retorna (a, b) desse intervalo
"""
def Plotar_e_encontrar_raiz(xmin=-2, xmax=2, espacamento=0.01):
    xs = []
    ys = []
    x = xmin
    while x <= xmax:
        xs.append(x)
        ys.append(F(x))
        x += espacamento
    plt.figure(figsize=(8,5))
    plt.axhline(0, color='black', linewidth=1)
    plt.plot(xs, ys, label="F(x)")
    plt.title("Gráfico de F(x) = sin(x²) + 1.1 - e^(-x)")
    plt.xlabel("x")
    plt.ylabel("F(x)")
    plt.legend()
    plt.grid(True)
    plt.show()
    for i in range(len(xs)-1):
        if ys[i] * ys[i+1] < 0:
            print(f"Raiz encontrada aproximadamente entre {xs[i]:.2f} e {xs[i+1]:.2f}")
            return xs[i], xs[i+1]
    print("Nenhuma raiz encontrada no intervalo visualizado.")
    sys.exit()

"""
Função principal, determina um subintervalo onde existe mudança de sinal, verifica se a raiz está contida 
num intervalo de largura 0.1, aplica ambos os métodos (Bisseção e Newton) e apresenta resultados obtidos.
Sem inputs
Retorna o subintervalo onde se encontra a raiz, aproximação por bisseção e Newton, erros majorados e números de iterações
"""
if __name__ == "__main__":
    print("=== Determinação automática do intervalo I ===")
    p, q = Plotar_e_encontrar_raiz(-2, 2)
    c = (p+q)/2
    a = c - 0.05
    b = c + 0.05
    print(f"Intervalo I definido automaticamente: [{a:.4f}, {b:.4f}] (amplitude = {b-a:.2f})")
    eps = 1e-9
    x0 = (a+b)/2
    print("\n Método da Bisseção ")
    x_bissecao, erro_b, it_b = Bissecao(a, b, eps)
    print(f"Resultado Bisseção: x = {x_bissecao:.12f} | erro <= {erro_b:.2e} | iterações = {it_b}")

    print("\nMétodo de Newton")
    x_newton, erro_n, it_n = Newton(a, b, eps, x0)
    print(f"Resultado Newton:   x = {x_newton:.12f} | erro <= {erro_n:.2e} | iterações = {it_n}")
