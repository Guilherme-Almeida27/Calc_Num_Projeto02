# Cálculo Numérico - Atividade Complementar (UNIVASF) - 2026.2
#
# Questões:
#   01.1 - Meta de produção ......... Eliminação de Gauss (método direto)
#   02.1 - Ponte de Wheatstone ...... Gauss-Seidel (método iterativo)
#   03.2 - Queda de voltagem ........ Interpolação de Lagrange e de Newton
#   04.2 - Área da superfície do rio  Trapézio e Simpson repetidos
#
# Como executar:  python atividade_calculo_numerico_v3.py
# Para os gráficos é preciso ter o matplotlib:  pip install matplotlib


try:
    import matplotlib.pyplot as plt
    tem_grafico = True
except ImportError:
    tem_grafico = False

LARGURA = 66        # largura das divisórias na tela


# ---------------------------------------------------------------
# FUNÇÕES DE FORMATAÇÃO (só cuidam da aparência da saída)
# ---------------------------------------------------------------

def titulo(texto):
    # caixa dupla em volta do título
    print("\n╔" + "═" * (LARGURA - 2) + "╗")
    print("║" + texto.center(LARGURA - 2) + "║")
    print("╚" + "═" * (LARGURA - 2) + "╝")


def secao(texto):
    # linha de divisão com o nome da seção:  ── Texto ─────────────
    print("\n── " + texto + " " + "─" * max(3, LARGURA - len(texto) - 4))


def aviso(texto):
    print(f"\n  [!] {texto}")


def tabela(cabecalhos, linhas, alinhar=None, recuo=2):
    # Desenha uma tabela com bordas. Todas as células já chegam como texto.
    # alinhar: um caractere por coluna  ("<" esquerda, ">" direita, "^" centro)
    n = len(cabecalhos)
    if alinhar is None:
        alinhar = ">" * n
    larguras = []
    for j in range(n):
        larguras.append(len(str(cabecalhos[j])))
    for linha in linhas:
        for j in range(n):
            larguras[j] = max(larguras[j], len(str(linha[j])))

    espaco = " " * recuo

    def divisoria(esq, meio, dir):
        return espaco + esq + meio.join("─" * (w + 2) for w in larguras) + dir

    def montar(celulas, alinhamentos):
        partes = []
        for j in range(n):
            partes.append(f" {str(celulas[j]):{alinhamentos[j]}{larguras[j]}} ")
        return espaco + "│" + "│".join(partes) + "│"

    print(divisoria("┌", "┬", "┐"))
    print(montar(cabecalhos, "^" * n))
    print(divisoria("├", "┼", "┤"))
    for linha in linhas:
        print(montar(linha, alinhar))
    print(divisoria("└", "┴", "┘"))


# ---------------------------------------------------------------
# FUNÇÕES AUXILIARES (leitura de dados)
# ---------------------------------------------------------------

def ler_numero(texto):
    # aceita vírgula ou ponto e repete a pergunta se o valor for inválido
    while True:
        try:
            return float(input(texto).replace(",", "."))
        except ValueError:
            print("    Valor inválido, tente de novo.")


def ler_inteiro(texto, minimo):
    while True:
        try:
            valor = int(input(texto))
            if valor >= minimo:
                return valor
        except ValueError:
            pass
        print(f"    Digite um número inteiro maior ou igual a {minimo}.")


def ler_lista(texto, quantidade):
    # lê vários números na mesma linha, separados por espaço
    # exemplo: 0,25 0,75 1,25
    while True:
        partes = input(texto).split()
        if len(partes) != quantidade:
            print(f"    Digite exatamente {quantidade} valores separados por espaço.")
            continue
        try:
            numeros = []
            for p in partes:
                numeros.append(float(p.replace(",", ".")))
            return numeros
        except ValueError:
            print("    Há algum valor inválido, tente de novo.")


def perguntar(texto):
    # devolve True para "s" ou Enter, e False para "n"
    resposta = input(texto + " (s/n) [Enter = s]: ").strip().lower()
    return resposta != "n"


def quer_grafico():
    if not tem_grafico:
        print("\n  (Instale o matplotlib para ver o gráfico: pip install matplotlib)")
        return False
    return perguntar("\n  Mostrar o gráfico?")


def mostrar_sistema(A, b, nomes):
    # mostra o sistema com as colunas alinhadas:
    #   (1)  + 4 x1   + 3 x2   + 2 x3   =  960
    n = len(A)
    textos = []         # textos[i][j] = "4 x1" (ou vazio, se o coeficiente é 0)
    largura = 0
    for i in range(n):
        linha = []
        for j in range(len(A[i])):
            if A[i][j] == 0:
                linha.append("")
            else:
                t = f"{abs(A[i][j]):g} {nomes[j]}"
                linha.append(t)
                largura = max(largura, len(t))
        textos.append(linha)

    largura_b = max(len(f"{valor:g}") for valor in b)
    for i in range(n):
        partes = []
        for j in range(len(A[i])):
            if textos[i][j] == "":
                partes.append(" " * (largura + 2))
            else:
                sinal = "+" if A[i][j] > 0 else "-"
                partes.append(f"{sinal} {textos[i][j]:>{largura}}")
        print(f"    ({i + 1})  " + "  ".join(partes) + f"  =  {b[i]:>{largura_b}g}")


def mostrar_matriz(M, b):
    # mostra a matriz aumentada [M | b]
    for i in range(len(M)):
        valores = "  ".join(f"{v:10.4f}" for v in M[i])
        print(f"    │ {valores}  ┆ {b[i]:10.4f} │")


# ---------------------------------------------------------------
# TÓPICO 01 - ELIMINAÇÃO DE GAUSS (MÉTODO DIRETO)
# ---------------------------------------------------------------

def eliminacao_gauss(A, b, mostrar_passos):
    # 1) Triangularização: zera os elementos abaixo da diagonal.
    #    Pivoteamento parcial: em cada coluna, usa como pivô o maior
    #    valor absoluto, para evitar divisão por zero e reduzir erros.
    # 2) Retrossubstituição: resolve de baixo para cima.
    n = len(A)
    M = []
    for linha in A:
        M.append(linha[:])      # cópia, para não alterar a matriz original
    v = b[:]

    if mostrar_passos:
        secao("Matriz aumentada inicial  [A | b]")
        mostrar_matriz(M, v)

    for k in range(n - 1):
        # procura a linha com o maior |valor| na coluna k
        p = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > abs(M[p][k]):
                p = i
        if abs(M[p][k]) < 1e-12:
            aviso("Pivô nulo: o sistema não tem solução única.")
            return None

        if mostrar_passos:
            secao(f"Etapa {k + 1}  (coluna {k + 1})")
        if p != k:
            M[k], M[p] = M[p], M[k]
            v[k], v[p] = v[p], v[k]
            if mostrar_passos:
                print(f"    Troca da linha {k + 1} com a linha {p + 1} (pivô parcial)")

        # zera os elementos abaixo do pivô
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]           # multiplicador
            if mostrar_passos:
                print(f"    m{i + 1}{k + 1} = {m:8.4f}   ->   L{i + 1} = L{i + 1} - ({m:.4f}) * L{k + 1}")
            for j in range(k, n):
                M[i][j] = M[i][j] - m * M[k][j]
            v[i] = v[i] - m * v[k]

        if mostrar_passos:
            print()
            mostrar_matriz(M, v)

    if abs(M[n - 1][n - 1]) < 1e-12:
        aviso("Pivô nulo: o sistema não tem solução única.")
        return None

    # retrossubstituição: começa pela última equação
    x = [0] * n
    for i in range(n - 1, -1, -1):
        soma = 0
        for j in range(i + 1, n):
            soma += M[i][j] * x[j]
        x[i] = (v[i] - soma) / M[i][i]
    return x


def questao_01_1():
    titulo("01.1 - PRODUÇÃO DE TRANSISTORES, RESISTORES E CHIPS")
    while True:
        print()
        if perguntar("  Usar os dados do enunciado (960, 510, 610)?"):
            produtos = ["Transistor", "Resistor", "Chip"]
            materiais = ["Cobre", "Zinco", "Vidro"]
            # cada linha é um material; cada coluna é um produto
            A = [[4, 3, 2],     # cobre usado por transistor, resistor e chip
                 [1, 3, 1],     # zinco
                 [2, 1, 3]]     # vidro
            b = [960, 510, 610]
        elif perguntar("  Manter a tabela e mudar só as quantidades disponíveis?"):
            produtos = ["Transistor", "Resistor", "Chip"]
            materiais = ["Cobre", "Zinco", "Vidro"]
            A = [[4, 3, 2], [1, 3, 1], [2, 1, 3]]
            b = []
            for nome in materiais:
                b.append(ler_numero(f"    {nome} disponível: "))
        else:
            n = ler_inteiro("  Quantos produtos (= quantos materiais)? ", 2)
            produtos = []
            for j in range(n):
                produtos.append(input(f"    Nome do produto {j + 1}: "))
            materiais = []
            A = []
            b = []
            for i in range(n):
                nome = input(f"\n    Nome do material {i + 1}: ")
                materiais.append(nome)
                A.append(ler_lista(f"    Quanto de {nome} cada produto usa ({n} valores): ", n))
                b.append(ler_numero(f"    Quantidade de {nome} disponível: "))

        # nomes curtos das incógnitas: x1, x2, x3...
        secao("Incógnitas")
        nomes = []
        for j in range(len(produtos)):
            nomes.append("x" + str(j + 1))
            print(f"    x{j + 1} = quantidade de {produtos[j]}")
        secao("Sistema de equações (uma para cada material)")
        mostrar_sistema(A, b, nomes)

        secao("Resolução por Eliminação de Gauss")
        x = eliminacao_gauss(A, b, True)
        if x is not None:
            secao("RESULTADO")
            linhas = []
            for j in range(len(produtos)):
                linhas.append([nomes[j], produtos[j], f"{x[j]:.4f}"])
            tabela(["", "Produto", "Quantidade"], linhas, "^<>")

            # conferência: material usado x material disponível
            secao("Conferência dos materiais")
            linhas = []
            for i in range(len(materiais)):
                usado = 0
                for j in range(len(produtos)):
                    usado += A[i][j] * x[j]
                linhas.append([materiais[i], f"{usado:.4f}", f"{b[i]:g}"])
            tabela(["Material", "Usado", "Disponível"], linhas, "<>>")

            # avisa se a produção não é possível na prática
            for valor in x:
                if valor < -1e-9:
                    aviso("Há quantidade negativa: essa meta de produção é impossível.")
                    break
                if abs(valor - round(valor)) > 1e-6:
                    aviso("Há quantidade não inteira: na prática é preciso arredondar.")
                    break

            if quer_grafico():
                grafico_producao(produtos, materiais, A, b, x)

        if not perguntar("\n  Novo cálculo nesta questão?"):
            break


def grafico_producao(produtos, materiais, A, b, x):
    plt.figure(figsize=(11, 4.5))

    # ---- gráfico 1: quantidade de cada produto ----
    plt.subplot(1, 2, 1)
    plt.bar(produtos, x, color="tab:blue")
    for j in range(len(x)):
        plt.text(j, x[j], f"{x[j]:.1f}", ha="center", va="bottom")
    plt.title("Produção da semana")
    plt.ylabel("Unidades")

    # ---- gráfico 2: material usado por cada produto (barras empilhadas) ----
    # mostra que a soma do que cada produto usa é igual ao disponível
    plt.subplot(1, 2, 2)
    base = [0] * len(materiais)
    for j in range(len(produtos)):
        usado = []
        for i in range(len(materiais)):
            usado.append(A[i][j] * x[j])
        plt.bar(materiais, usado, bottom=base, label=produtos[j])
        for i in range(len(materiais)):
            base[i] += usado[i]
    plt.plot(materiais, b, "k_", markersize=40, markeredgewidth=2, label="Disponível")
    plt.title("Uso de cada material")
    plt.ylabel("Unidades de material")
    plt.ylim(0, max(b) * 1.4)       # deixa espaço para a legenda
    plt.legend(loc="upper right")

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------
# TÓPICO 02 - GAUSS-SEIDEL (MÉTODO ITERATIVO)
# ---------------------------------------------------------------

def sassenfeld(A):
    # Critério de Sassenfeld: se o maior beta for menor que 1,
    # o Gauss-Seidel converge com certeza.
    n = len(A)
    beta = [0] * n
    for i in range(n):
        soma = 0
        for j in range(n):
            if j < i:
                soma += abs(A[i][j]) * beta[j]   # usa os betas já calculados
            elif j > i:
                soma += abs(A[i][j])
        beta[i] = soma / abs(A[i][i])
    return max(beta)


def gauss_seidel(A, b, nomes, eps, max_iter):
    # x_i = (b_i - soma dos outros termos da linha i) / a_ii
    # usando sempre o valor mais novo de cada variável.
    # Para quando o erro relativo  max|x(k) - x(k-1)| / max|x(k)|  < eps
    n = len(A)
    x = [0] * n                     # aproximação inicial: vetor nulo
    historico = [x[:]]              # guarda x de cada iteração (para o gráfico)
    erros = []
    linhas = []                     # linhas da tabela de iterações
    convergiu = False

    for k in range(1, max_iter + 1):
        x_anterior = x[:]
        for i in range(n):
            soma = 0
            for j in range(n):
                if j != i:
                    soma += A[i][j] * x[j]
            x[i] = (b[i] - soma) / A[i][i]

        maior_diferenca = 0
        maior_valor = 0
        for i in range(n):
            maior_diferenca = max(maior_diferenca, abs(x[i] - x_anterior[i]))
            maior_valor = max(maior_valor, abs(x[i]))
        if maior_valor == 0:        # solução nula: evita divisão por zero
            erro = 0
        else:
            erro = maior_diferenca / maior_valor

        historico.append(x[:])
        erros.append(erro)
        linhas.append([str(k)] + [f"{valor:.6f}" for valor in x] + [f"{erro:.3e}"])
        if erro < eps:
            convergiu = True
            break

    # mostra a tabela (se for muito longa, mostra o começo e o fim)
    if len(linhas) > 30:
        linhas = linhas[:10] + [["..."] * (n + 2)] + linhas[-5:]
    tabela(["k"] + nomes + ["erro relativo"], linhas)

    if convergiu:
        print(f"\n  Convergiu em {len(erros)} iterações (precisão {eps:g}).")
    else:
        aviso("Não convergiu no número máximo de iterações.")
    return x, historico, erros


def questao_02_1():
    titulo("02.1 - PONTE DE WHEATSTONE (GAUSS-SEIDEL)")
    while True:
        print()
        if perguntar("  Usar os dados do enunciado (E=30 V, R1=20, R2=R3=R4=R5=120)?"):
            E, R1, R2, R3, R4, R5 = 30, 20, 120, 120, 120, 120
        else:
            E = ler_numero("    E (Volts) = ")
            R1 = ler_numero("    R1 (Ohms) = ")
            R2 = ler_numero("    R2 (Ohms) = ")
            R3 = ler_numero("    R3 (Ohms) = ")
            R4 = ler_numero("    R4 (Ohms) = ")
            R5 = ler_numero("    R5 (Ohms) = ")
            if min(R1, R2, R3, R4, R5) <= 0:
                aviso("As resistências devem ser maiores que zero.")
                continue
        eps = ler_numero("    Precisão (ex: 0,0001) = ")
        if eps <= 0:
            aviso("A precisão deve ser maior que zero.")
            continue

        secao("Leis de Kirchhoff (sentidos das correntes como na figura)")
        print("    Nó A:            I1 + I2 = I6")
        print("    Nó B:            I1 = I4 + I5")
        print("    Nó C:            I2 + I5 = I3")
        print("    Malha ABC:       R1*I1 + R5*I5 - R2*I2 = 0")
        print("    Malha BCD:       R5*I5 + R3*I3 - R4*I4 = 0")
        print("    Malha da fonte:  R1*I1 + R4*I4 = E")

        # Essas 6 equações têm zeros na diagonal e nenhuma ordem delas
        # garante que o Gauss-Seidel converge (com R1=R2=R3=R4=R5 ele não converge em ordem nenhuma). Por isso usamos correntes de malha:
        #   J1 = malha A-C-B,  J2 = malha B-C-D,  J3 = malha da fonte
        # Cada corrente da figura é uma combinação delas:
        #   I1 = J3 - J1   I2 = J1   I3 = J2   I4 = J3 - J2   I5 = J2 - J1   I6 = J3
        # Substituindo nas equações das malhas, sobra um sistema 3x3 em que a
        #diagonal "domina": o Gauss-Seidel converge para quaisquer R > 0.

        A = [[R1 + R2 + R5, -R5, -R1],
             [-R5, R3 + R4 + R5, -R4],
             [-R1, -R4, R1 + R4]]
        b = [0, 0, E]

        secao("Sistema com as correntes de malha J1, J2, J3")
        print("    I1 = J3 - J1   I2 = J1   I3 = J2   I4 = J3 - J2   I5 = J2 - J1   I6 = J3\n")
        mostrar_sistema(A, b, ["J1", "J2", "J3"])

        beta = sassenfeld(A)
        if beta < 1:
            print(f"\n  Critério de Sassenfeld: maior beta = {beta:.4f} (< 1, converge)")
        else:
            aviso(f"Critério de Sassenfeld: maior beta = {beta:.4f} (>= 1, convergência não garantida)")
        
        secao("Iterações do Gauss-Seidel (começando com J = 0)")
        J, historico, erros = gauss_seidel(A, b, ["J1", "J2", "J3"], eps, 1000)

        # volta para as correntes da figura
        correntes = de_malha_para_ramo(J)

        # conferência: resolve as 6 equações originais pelo método direto
        A6 = [[1,  1, 0, 0,   0, -1],
              [1,  0, 0, -1,  -1, 0],
              [ 0,  1, -1, 0,  1, 0],
              [R1, -R2, 0, 0, R5, 0],
              [0, 0, R3, -R4, R5, 0],
              [R1,  0, 0,  R4, 0, 0]]
        
        b6 = [0, 0, 0, 0, 0, E]
        exata = eliminacao_gauss(A6, b6, False)

        secao("CORRENTES (Gauss-Seidel x conferência pelo método direto)")
        linhas = []
        for i in range(6):
            if exata is None:
                linhas.append([f"I{i + 1}", f"{correntes[i]:.6f}", f"{correntes[i] * 1000:.3f}", "-", "-"])
            else:
                linhas.append([f"I{i + 1}", f"{correntes[i]:.6f}", f"{correntes[i] * 1000:.3f}",
                               f"{exata[i]:.6f}", f"{abs(exata[i] - correntes[i]):.2e}"])
        tabela(["", "G-Seidel (A)", "G-Seidel (mA)", "Direto (A)", "Diferença"], linhas, "^>>>>")

        if quer_grafico():
            grafico_wheatstone(historico, erros)

        if not perguntar("\n  Novo cálculo nesta questão?"):
            break


def de_malha_para_ramo(J):
    J1, J2, J3 = J
    return [J3 - J1, J1, J2, J3 - J2, J2 - J1, J3]


def grafico_wheatstone(historico, erros):
    plt.figure(figsize=(11, 4.5))

    # ---- gráfico 1: valor de cada corrente a cada iteração ----
    plt.subplot(1, 2, 1)
    iteracoes = range(len(historico))
    for i in range(6):
        valores = []
        for J in historico:
            valores.append(de_malha_para_ramo(J)[i])
        plt.plot(iteracoes, valores, marker=".", label=f"I{i + 1}")
    plt.title("Correntes a cada iteração")
    plt.xlabel("Iteração")
    plt.ylabel("Corrente (A)")
    plt.grid(True)
    plt.legend()

    # ---- gráfico 2: erro relativo a cada iteração ----
    plt.subplot(1, 2, 2)
    plt.plot(range(1, len(erros) + 1), erros, "o-")
    plt.yscale("log")       # escala logarítmica, porque o erro cai muito rápido
    plt.title("Velocidade de convergência")
    plt.xlabel("Iteração")
    plt.ylabel("Erro relativo")
    plt.grid(True)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------
# TÓPICO 03 - INTERPOLAÇÃO (LAGRANGE E NEWTON)
# ---------------------------------------------------------------

def escolher_pontos(xs, ys, x_procurado, quantidade):
    # pega os pontos mais próximos do valor procurado
    # (o enunciado pede pontos-base em torno do valor desconhecido)
    pares = []
    for i in range(len(xs)):
        distancia = abs(xs[i] - x_procurado)
        pares.append((distancia, xs[i], ys[i]))
    pares.sort()                    # ordena pela distância
    escolhidos = pares[:quantidade]

    # coloca os escolhidos em ordem crescente de x
    pontos = []
    for par in escolhidos:
        pontos.append((par[1], par[2]))     # (x, y)
    pontos.sort()                           # ordena pelo x
    px = []
    py = []
    for ponto in pontos:
        px.append(ponto[0])
        py.append(ponto[1])
    return px, py


def lagrange(xs, ys, x, mostrar_passos):
    # P(x) = soma de y_k * L_k(x)
    # L_k(x) = produto de (x - x_j) / (x_k - x_j), para todo j diferente de k
    resultado = 0
    linhas = []
    for k in range(len(xs)):
        L = 1
        for j in range(len(xs)):
            if j != k:
                L = L * (x - xs[j]) / (xs[k] - xs[j])
        linhas.append([str(k), f"{xs[k]:g}", f"{ys[k]:g}", f"{L:.6f}", f"{ys[k] * L:.6f}"])
        resultado += ys[k] * L
    if mostrar_passos:
        tabela(["k", "xk", "yk", f"Lk({x:g})", "yk * Lk"], linhas, "^>>>>", 4)
    return resultado


def newton(xs, ys, x, mostrar_passos):
    # Tabela de diferenças divididas:
    # tabela[0] são os y, tabela[1] as diferenças de ordem 1, e assim por diante.
    n = len(xs)
    dif = [ys[:]]
    for ordem in range(1, n):
        anterior = dif[ordem - 1]
        nova = []
        for i in range(n - ordem):
            nova.append((anterior[i + 1] - anterior[i]) / (xs[i + ordem] - xs[i]))
        dif.append(nova)

    if mostrar_passos:
        linhas = []
        for i in range(n):
            linha = [f"{xs[i]:g}"]
            for ordem in range(n):
                if i < n - ordem:
                    linha.append(f"{dif[ordem][i]:.6f}")
                else:
                    linha.append("")
            linhas.append(linha)
        tabela(["xi"] + [f"ordem {o}" for o in range(n)], linhas, ">" * (n + 1), 4)
        print("    Coeficientes (1ª linha de cada ordem):")
        for inicio in range(0, n, 3):       # no máximo 3 por linha
            coef = "   ".join(f"d{k} = {dif[k][0]:10.6f}" for k in range(inicio, min(inicio + 3, n)))
            print("      " + coef)

    # P(x) = d0 + d1*(x-x0) + d2*(x-x0)*(x-x1) + ...
    resultado = 0
    produto = 1
    for k in range(n):
        resultado += dif[k][0] * produto
        produto = produto * (x - xs[k])
    return resultado


def questao_03_2():
    titulo("03.2 - INTERPOLAÇÃO (LAGRANGE E NEWTON)")
    while True:
        print()
        if perguntar("  Usar os dados do enunciado (estimar V para i = 1,15)?"):
            xs = [0.25, 0.75, 1.25, 1.5, 2.0]
            ys = [-0.45, -0.60, 0.70, 1.88, 6.0]
            x = 1.15
        else:
            n = ler_inteiro("    Quantos pontos tem a tabela (mínimo 3)? ", 3)
            xs = ler_lista("    Valores de x: ", n)
            if len(set(xs)) != n:       # set() remove os repetidos
                aviso("Os valores de x não podem se repetir.")
                continue
            ys = ler_lista("    Valores de y: ", n)
            x = ler_numero("    Valor de x para estimar y: ")

        secao("Dados")
        linhas = []
        for i in range(len(xs)):
            linhas.append([f"{xs[i]:g}", f"{ys[i]:g}"])
        tabela(["x", "y"], linhas)
        print(f"\n  Valor a estimar: x = {x:g}")

        if x < min(xs) or x > max(xs):
            aviso("O valor está fora da tabela (extrapolação).")

        grau_maximo = min(4, len(xs) - 1)
        resultados = []
        resumo = []
        for grau in range(2, grau_maximo + 1):
            px, py = escolher_pontos(xs, ys, x, grau + 1)
            pontos_texto = "  ".join(f"{p:g}" for p in px)
            secao(f"GRAU {grau}  -  pontos-base: {pontos_texto}")
            print("\n  Forma de Lagrange:")
            valor_lagrange = lagrange(px, py, x, True)
            print(f"    P{grau}({x:g}) = soma de yk * Lk = {valor_lagrange:.6f}")
            print("\n  Forma de Newton:")
            valor_newton = newton(px, py, x, True)
            print(f"    P{grau}({x:g}) = {valor_newton:.6f}")
            resultados.append([grau, valor_newton, px, py])
            resumo.append([str(grau), pontos_texto, f"{valor_lagrange:.6f}", f"{valor_newton:.6f}"])

        secao("RESUMO")
        tabela(["Grau", "Pontos-base", "Lagrange", "Newton"], resumo, "^<>>")
        print("\n  Lagrange e Newton dão o mesmo valor: é o mesmo polinômio")
        print("  escrito de formas diferentes.")

        if quer_grafico():
            grafico_interpolacao(xs, ys, x, resultados)

        if not perguntar("\n  Novo cálculo nesta questão?"):
            break


def grafico_interpolacao(xs, ys, x, resultados):
    plt.figure(figsize=(8, 5))
    plt.plot(xs, ys, "ko", markersize=8, label="Dados da tabela")

    # desenha cada polinômio com 200 pontos entre o menor e o maior x
    inicio = min(min(xs), x)
    fim = max(max(xs), x)
    for grau, valor, px, py in resultados:
        pontos_x = []
        pontos_y = []
        for i in range(201):
            t = inicio + i * (fim - inicio) / 200
            pontos_x.append(t)
            pontos_y.append(lagrange(px, py, t, False))
        plt.plot(pontos_x, pontos_y, label=f"Grau {grau}: P({x:g}) = {valor:.4f}")
        plt.plot(x, valor, "x", markersize=10, color="red")

    plt.axvline(x, color="gray", linestyle=":")
    plt.title("Interpolação polinomial")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()


# ---------------------------------------------------------------
# TÓPICO 04 - INTEGRAÇÃO NUMÉRICA (TRAPÉZIO E SIMPSON REPETIDOS)
# ---------------------------------------------------------------

def trapezio_repetido(y, h):
    # I = h/2 * [y0 + 2*(y1 + y2 + ... + y(n-1)) + yn]
    n = len(y) - 1
    meio = 0
    for i in range(1, n):
        meio += y[i]
    return h / 2 * (y[0] + 2 * meio + y[n])


def simpson_repetido(y, h):
    # I = h/3 * [y0 + 4*(soma dos ímpares) + 2*(soma dos pares) + yn]
    # Só funciona com número PAR de subintervalos.
    n = len(y) - 1
    if n % 2 != 0:
        return None
    impares = 0
    pares = 0
    for i in range(1, n):
        if i % 2 == 1:
            impares += y[i]
        else:
            pares += y[i]
    return h / 3 * (y[0] + 4 * impares + 2 * pares + y[n])


def questao_04_2():
    titulo("04.2 - ÁREA DA SUPERFÍCIE DO RIO")
    print("\n  Área = integral de [M2(x) - M1(x)] dx")
    while True:
        print()
        if perguntar("  Usar os dados do enunciado?"):
            xs = [0, 10, 20, 30, 40]
            m1 = [50.8, 86.2, 136, 72.8, 51]
            m2 = [113.6, 144.5, 185, 171.2, 95.3]
        else:
            n = ler_inteiro("    Quantos pontos (mínimo 3)? ", 3)
            xs = ler_lista("    Valores de x (crescentes e igualmente espaçados): ", n)
            m1 = ler_lista("    Valores de M1(x): ", n)
            m2 = ler_lista("    Valores de M2(x): ", n)

        # confere se os pontos são crescentes e igualmente espaçados
        h = xs[1] - xs[0]
        espacamento_ok = h > 0
        for i in range(len(xs) - 1):
            if abs((xs[i + 1] - xs[i]) - h) > 1e-9:
                espacamento_ok = False
        if not espacamento_ok:
            aviso("Os valores de x devem ser crescentes e igualmente espaçados.")
            continue

        # largura do rio em cada ponto = distância entre as margens
        largura = []
        for i in range(len(xs)):
            largura.append(m2[i] - m1[i])

        secao("Dados")
        linhas = []
        for i in range(len(xs)):
            linhas.append([f"{xs[i]:g}", f"{m1[i]:.2f}", f"{m2[i]:.2f}", f"{largura[i]:.2f}"])
        tabela(["x (m)", "M1(x)", "M2(x)", "Largura"], linhas)
        print(f"\n  h = {h:g} m      n = {len(xs) - 1} subintervalos")

        area_trapezio = trapezio_repetido(largura, h)
        area_simpson = simpson_repetido(largura, h)

        secao("RESULTADO")
        linhas = [["Trapézio repetido", f"{area_trapezio:.4f}"]]
        if area_simpson is None:
            linhas.append(["Simpson repetido", "não se aplica"])
        else:
            linhas.append(["Simpson repetido", f"{area_simpson:.4f}"])
            linhas.append(["Diferença", f"{abs(area_simpson - area_trapezio):.4f}"])
        tabela(["Método", "Área (m²)"], linhas, "<>")
        if area_simpson is None:
            print("\n  Simpson exige número par de subintervalos (número ímpar aqui).")

        if quer_grafico():
            grafico_rio(xs, m1, m2, largura, area_trapezio, area_simpson)

        if not perguntar("\n  Novo cálculo nesta questão?"):
            break


def grafico_rio(xs, m1, m2, largura, area_trapezio, area_simpson):
    plt.figure(figsize=(11, 4.5))

    # ---- gráfico 1: as margens do rio ----
    plt.subplot(1, 2, 1)
    plt.plot(xs, m1, "o-", label="Margem M1")
    plt.plot(xs, m2, "s-", label="Margem M2")
    plt.fill_between(xs, m1, m2, alpha=0.3, label="Rio")
    plt.title("Superfície do rio")
    plt.xlabel("x (m)")
    plt.ylabel("Distância à linha de referência (m)")
    plt.grid(True)
    plt.legend()

    # ---- gráfico 2: a largura e os trapézios usados na regra ----
    plt.subplot(1, 2, 2)
    plt.fill_between(xs, largura, alpha=0.3, color="tab:orange",
                     label=f"Trapézios: {area_trapezio:.1f} m²")
    plt.plot(xs, largura, "o-", color="tab:orange")
    for i in range(len(xs)):
        plt.plot([xs[i], xs[i]], [0, largura[i]], color="tab:orange", linewidth=0.8)
    if area_simpson is not None:
        # parábolas da regra de Simpson: uma a cada 2 subintervalos
        for i in range(0, len(xs) - 2, 2):
            px = xs[i:i + 3]
            py = largura[i:i + 3]
            curva_x = []
            curva_y = []
            for k in range(51):
                t = px[0] + k * (px[2] - px[0]) / 50
                curva_x.append(t)
                curva_y.append(lagrange(px, py, t, False))
            rotulo = f"Simpson: {area_simpson:.1f} m²" if i == 0 else None
            plt.plot(curva_x, curva_y, "--", color="tab:green", label=rotulo)
    plt.title("Largura do rio (M2 - M1)")
    plt.xlabel("x (m)")
    plt.ylabel("Largura (m)")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------
# PROGRAMA PRINCIPAL
# ---------------------------------------------------------------

def mostrar_menu():
    print("\n╔" + "═" * (LARGURA - 2) + "╗")
    print("║" + "CÁLCULO NUMÉRICO - ATIVIDADE COMPLEMENTAR".center(LARGURA - 2) + "║")
    print("╠" + "═" * (LARGURA - 2) + "╣")
    itens = ["1  -  01.1  Meta de produção          (Eliminação de Gauss)",
             "2  -  02.1  Ponte de Wheatstone       (Gauss-Seidel)",
             "3  -  03.2  Queda de voltagem         (Lagrange e Newton)",
             "4  -  04.2  Área do rio               (Trapézio e Simpson)",
             "0  -  Sair"]
    for item in itens:
        print("║  " + item.ljust(LARGURA - 4) + "║")
    print("╚" + "═" * (LARGURA - 2) + "╝")


if __name__ == "__main__":

    print("\n  Números decimais podem ser digitados com vírgula ou ponto.")
    print("  Em listas, separe os valores por espaço.  Ex: 0,25 0,75 1,25")

    opcao = ""
    while opcao != "0":
        mostrar_menu()
        opcao = input("  Escolha uma opção: ").strip()

        # o try evita que um erro inesperado feche o programa
        try:
            if opcao == "1":
                questao_01_1()
            elif opcao == "2":
                questao_02_1()
            elif opcao == "3":
                questao_03_2()
            elif opcao == "4":
                questao_04_2()
            elif opcao != "0":
                aviso("Opção inválida.")
        except KeyboardInterrupt:
            print("\n  Cálculo interrompido.")
        except Exception as erro:
            print("\n  Erro:", erro)
            print("  Confira os valores digitados.")

    print("\n  Programa encerrado.\n")
