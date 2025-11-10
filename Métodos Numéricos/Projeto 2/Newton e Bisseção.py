import math
import sys

"""
Documentação deste ficheiro, para melhor compreensão, foi feita de acordo com a seguinte estrutura:
Explicação da função, Inputs, Outputs
Comentários removidos após última análise e a documentação aparece antes da função
"""

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
    if fa*fb >= 0:
        raise ValueError("F(a) e F(b) devem ter sinais opostos.")   
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
Recebe os valores de a e b (float), tem predefenido 2000 samples (int), que funciona para dividir o intervalo dado em 2000 partes
Retorna o valor do majorante do erro (float)
"""
def Majorante_erro(a, b, samples = 2000):
    espacamento = (b-a)/samples
    m_erro = float('inf')
    x = a
    for _ in range(samples + 1):
        m_erro = min(m_erro, abs(F_derivada(x)))
        x += espacamento
    if m_erro < 1*(10**(-14)):
        m_erro = 1*(10**(-14))
    return m_erro

"""
Aplicação do método de Newton
Recebe os valores de a, b, eps, todos float
Retorna os valor aproximado da raiz (float), o valor do erro majorado (float), número de iterações (int)
"""
def Newton(a, b, eps, x0=None):
    if x0 is None:
        x0 = (a+b)/2
    m_erro = Majorante_erro(a, b)
    x = x0
    max_iteracoes = 5
    iteracoes = 0    
    while iteracoes < max_iteracoes:
        fx = F(x)
        dfx = F_derivada(x)
        if abs(dfx) < (1*(10**(-14))):
            x1 = (a+b)/2
        else:
            x1 = x-(fx/dfx)
        if not (a <= x1 <= b):
            x1 = (a+b)/2
        iteracoes += 1
        erro_majorado = abs(F(x1)) / m_erro
        if erro_majorado <= eps:
            return x1, erro_majorado, iteracoes
        x = x1
    return x, abs(F(x))/m_erro, iteracoes

"""
Verifica a existência e localização de uma raiz no intervalo dado, iterando pela incrementação de k
Recebe um intervalo, xmin e xmax, ambos float
Retorna um subintervalo onde occore a mudança de sinal, ambos float, se ocorrer, senão retorna None
"""
def Raiz(xmin, xmax):
    k = 0.01
    while (xmin+k) <= xmax:
        a = F(xmin)
        b = F(xmin+k)
        if (a*b) < 0:
            print("Raiz encontrada para o intervalo:", xmin, xmax)
            return xmin, xmin + k
        else:
            xmin += k
    print("Não foi encontrada raiz no intervalo fornecido.")
    return None

"""
Função principal, lê o intervalo dado pelo utilizador, determina um subintervalo onde existe mudança de sinal,
verifica se a raiz está contida num intervalo de largura 0.1, aplica ambos os métodos (Bisseção e Newton) e apresenta resultados.
Obtém os valores xmin e xmax através do terminal, ambos float
Retorna o subintervalo onde se encontra a raiz, aproximação por bisseção e Newton, erros majorados e números de iterações
"""
if __name__ == "__main__":
    xmin = float(input("xmin = "))
    xmax = float(input("xmax = "))
    resultado = Raiz(xmin, xmax)
    if resultado is None:
        sys.exit("Não foi encontrada raiz no intervalo fornecido.")
    p, q = resultado
    c = (p+q)/2
    a = c - 0.05
    b = c + 0.05
    if (a <= p) and (q <= b):
        print("Verificação: o intervalo contém a raiz.")
    else:
        sys.exit("Verificação: o intervalo não contém a raiz, processo terminado.")
    eps = 1*(10**(-9))
    x0 = (a+b)/2
    Fa = F(a)
    Fb = F(b)
    if Fa*Fb < 0:
        a_bissecao, b_bissecao = a, b
        print("Ocorre mudança de sinal, valores para bisseção serão:", a, b)
    else:
        a_bissecao, b_bissecao = p, q
        print("Não ocorre mudança de sinal, valores para bisseção serão:", p, q)
    x_bissecao, erro_majorado_bissecao, iteracoes_bissecao = Bissecao(a_bissecao, b_bissecao, eps)
    print(f"Bisseção | x = {x_bissecao:.12f} | erro <= {erro_majorado_bissecao:.2e} | iterações = {iteracoes_bissecao}")
    print("Os valores usados para o método de Newton serão:", a, b, eps, x0)
    x_newton, erro_majorado_newton, iteracoes_newton = Newton(a, b, eps, x0)
    print(f"Newton | x = {x_newton:.12f} | erro <= {erro_majorado_newton:.2e} | iterações = {iteracoes_newton}")
