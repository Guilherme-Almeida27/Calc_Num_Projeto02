# Cálculo Numérico – Atividade Complementar 2026.2

Programa em Python que resolve quatro problemas de engenharia com métodos numéricos, desenvolvido para a disciplina de **Cálculo Numérico** da **UNIVASF**.

Além de calcular, o programa mostra os passos de cada método (etapas da eliminação, tabela de iterações, diferenças divididas…), aceita novos dados para resolver problemas semelhantes e oferece gráficos para facilitar a interpretação dos resultados.

## Problemas resolvidos

| Questão | Problema | Método numérico |
|:---:|---|---|
| 01.1 | Meta de produção de transistores, resistores e chips | Eliminação de Gauss com pivoteamento parcial |
| 02.1 | Correntes em uma Ponte de Wheatstone | Gauss-Seidel |
| 03.2 | Queda de voltagem em um resistor | Interpolação de Lagrange e de Newton |
| 04.2 | Área da superfície de um rio entre duas margens | Regras do Trapézio e de Simpson repetidas |

## Funcionalidades

- **Menu interativo**: escolha a questão, faça quantos cálculos quiser e saia quando terminar.
- **Dados do enunciado ou novos dados**: cada questão pergunta se deve usar os valores do enunciado (basta apertar Enter) ou pedir valores novos.
- **Passo a passo na tela**: matrizes de cada etapa, tabelas de iterações, tabelas de Lagrange e de diferenças divididas.
- **Conferência dos resultados**: o Gauss-Seidel é conferido pelo método direto e a produção é conferida pelo material consumido.
- **Gráficos** (opcionais, com matplotlib) em todas as questões.
- **Validação da entrada**: aceita vírgula ou ponto como separador decimal, repete a pergunta quando o valor é inválido e avisa sobre dados impossíveis (x repetidos, espaçamento desigual, resistência ≤ 0…).
- **Código simples**: só `for`, `if` e listas, sem bibliotecas de álgebra linear; todos os métodos foram implementados à mão, com comentários explicando cada fórmula.

## Requisitos

- Python 3.8 ou superior (testado no Python 3.13)
- [matplotlib](https://matplotlib.org/) — **opcional**, só para os gráficos. Sem ele, o programa funciona normalmente e apenas não mostra os gráficos.

## Como executar

```bash
# 1. Baixe o projeto
git clone https://github.com/SEU-USUARIO/SEU-REPOSITORIO.git
cd SEU-REPOSITORIO

# 2. (Opcional) instale o matplotlib para ver os gráficos
pip install matplotlib

# 3. Execute
python atividade_calculo_numerico_v3.py
```

> No Linux/macOS, use `python3` no lugar de `python` se for o caso.

## Como usar

Ao iniciar, o programa mostra o menu principal:

```
╔════════════════════════════════════════════════════════════════╗
║           CÁLCULO NUMÉRICO - ATIVIDADE COMPLEMENTAR            ║
╠════════════════════════════════════════════════════════════════╣
║  1  -  01.1  Meta de produção          (Eliminação de Gauss)   ║
║  2  -  02.1  Ponte de Wheatstone       (Gauss-Seidel)          ║
║  3  -  03.2  Queda de voltagem         (Lagrange e Newton)     ║
║  4  -  04.2  Área do rio               (Trapézio e Simpson)    ║
║  0  -  Sair                                                    ║
╚════════════════════════════════════════════════════════════════╝
```

Dicas para digitar os dados:

- Números decimais podem usar **vírgula ou ponto**: `0,25` ou `0.25`.
- Listas são digitadas **na mesma linha, separadas por espaço**: `0,25 0,75 1,25 1,5 2`.
- Nas perguntas de sim/não, **Enter equivale a "sim"**.

---

## Os problemas em detalhe

### 01.1 – Meta de produção (Eliminação de Gauss)

Uma fábrica produz transistores, resistores e chips usando cobre, zinco e vidro. Cada material vira uma equação (o total consumido deve ser igual ao disponível) e cada produto vira uma incógnita:

```
4·x1 + 3·x2 + 2·x3 = 960    (cobre)
1·x1 + 3·x2 + 1·x3 = 510    (zinco)
2·x1 + 1·x2 + 3·x3 = 610    (vidro)
```

O sistema é resolvido por **eliminação de Gauss com pivoteamento parcial**, mostrando a matriz aumentada, os multiplicadores e a matriz após cada etapa.

**Resultado:** 120 transistores, 100 resistores e 90 chips.

```
  ┌────┬────────────┬────────────┐
  │    │  Produto   │ Quantidade │
  ├────┼────────────┼────────────┤
  │ x1 │ Transistor │   120.0000 │
  │ x2 │ Resistor   │   100.0000 │
  │ x3 │ Chip       │    90.0000 │
  └────┴────────────┴────────────┘
```

Também é possível mudar só as quantidades da semana ou montar uma tabela nova, com qualquer número de produtos. O programa avisa quando a solução tem quantidade negativa (meta impossível) ou não inteira.

**Gráfico:** produção de cada item e uso de cada material (barras empilhadas), comparado com a quantidade disponível.

### 02.1 – Ponte de Wheatstone (Gauss-Seidel)

As Leis de Kirchhoff geram 6 equações para as correntes I1 a I6. Esse sistema tem zeros na diagonal, e o Gauss-Seidel não tem convergência garantida nele. Por isso o programa usa **correntes de malha** (J1, J2, J3), que reduzem o problema a um sistema 3×3 diagonalmente dominante:

```
(R1+R2+R5)·J1 −        R5·J2 −        R1·J3 = 0
      −R5·J1 + (R3+R4+R5)·J2 −        R4·J3 = 0
      −R1·J1 −        R4·J2 +   (R1+R4)·J3 = E
```

As correntes da figura são obtidas por: I1 = J3 − J1, I2 = J1, I3 = J2, I4 = J3 − J2, I5 = J2 − J1 e I6 = J3.

O programa verifica o **critério de Sassenfeld**, mostra a tabela de iterações com o erro relativo e confere o resultado resolvendo as 6 equações originais por eliminação de Gauss.

**Resultado** (E = 30 V, R1 = 20 Ω, R2 = R3 = R4 = R5 = 120 Ω, precisão 0,0001): convergência em 14 iterações, com β de Sassenfeld = 0,5385.

| Corrente | Valor (A) |
|:---:|:---:|
| I1 | 0,2609 |
| I2 | 0,0978 |
| I3 | 0,1521 |
| I4 | 0,2065 |
| I5 | 0,0543 |
| I6 | 0,3587 |

**Gráfico:** valor de cada corrente a cada iteração e queda do erro (escala logarítmica).

### 03.2 – Queda de voltagem (Lagrange e Newton)

A partir da tabela de corrente (i) e voltagem (V) medidas em um resistor, o programa estima V para **i = 1,15** com polinômios de grau 2, 3 e 4. Os pontos-base são escolhidos automaticamente entre os mais próximos do valor procurado.

Para cada grau, são mostradas a tabela dos L_k(x) de Lagrange e a tabela de diferenças divididas de Newton.

**Resultado:**

| Grau | Pontos-base | V(1,15) |
|:---:|---|:---:|
| 2 | 0,75 · 1,25 · 1,5 | 0,326933 |
| 3 | 0,75 · 1,25 · 1,5 · 2 | 0,347840 |
| 4 | 0,25 · 0,75 · 1,25 · 1,5 · 2 | 0,337286 |

Lagrange e Newton dão o mesmo valor em cada grau, porque geram o mesmo polinômio escrito de formas diferentes. O programa avisa quando o valor está fora da tabela (extrapolação).

**Gráfico:** pontos da tabela, os três polinômios e as estimativas.

### 04.2 – Área da superfície do rio (Trapézio e Simpson)

A largura do rio em cada ponto é a distância entre as margens, M2(x) − M1(x). A área é a integral dessa largura no trecho [0, 40] m, calculada pelas **regras do Trapézio e de Simpson repetidas** (h = 10 m, n = 4 subintervalos).

**Resultado:**

| Método | Área |
|---|:---:|
| Trapézio repetido | 2592,5 m² |
| Simpson repetido | 2773,0 m² |

A regra de Simpson aproxima a largura por parábolas e acompanha melhor o aumento da largura em x = 30 m, por isso dá uma área maior. O programa confere se os pontos estão em ordem e igualmente espaçados, e avisa quando o número de subintervalos é ímpar (Simpson não se aplica).

**Gráfico:** as margens com a superfície do rio preenchida, e a largura com os trapézios e as parábolas de Simpson desenhados.

---

## Estrutura do código

Tudo fica em um único arquivo, `atividade_calculo_numerico_v3.py`, dividido em seções:

| Seção | Funções principais | O que faz |
|---|---|---|
| Formatação | `titulo`, `secao`, `aviso`, `tabela`, `mostrar_sistema`, `mostrar_matriz` | Aparência da saída (caixas, tabelas com bordas, sistemas alinhados) |
| Leitura de dados | `ler_numero`, `ler_inteiro`, `ler_lista`, `perguntar` | Entrada com validação |
| Tópico 01 | `eliminacao_gauss`, `questao_01_1`, `grafico_producao` | Método direto |
| Tópico 02 | `sassenfeld`, `gauss_seidel`, `questao_02_1`, `de_malha_para_ramo`, `grafico_wheatstone` | Método iterativo |
| Tópico 03 | `escolher_pontos`, `lagrange`, `newton`, `questao_03_2`, `grafico_interpolacao` | Interpolação |
| Tópico 04 | `trapezio_repetido`, `simpson_repetido`, `questao_04_2`, `grafico_rio` | Integração numérica |
| Programa principal | `mostrar_menu` | Menu e laço principal |

Os métodos numéricos ficam separados da entrada e da saída. Por exemplo, `eliminacao_gauss` é usada tanto na questão 01.1 quanto na conferência da questão 02.1, e `lagrange` é reaproveitada para desenhar as parábolas de Simpson no gráfico da questão 04.2.

## Limitações conhecidas

- **Gauss-Seidel com resistências muito diferentes:** quando os valores das resistências diferem em várias ordens de grandeza (por exemplo, R5 = 10 kΩ mantendo os demais valores do enunciado), o critério de Sassenfeld fica próximo de 1. A convergência fica lenta e o critério de parada pode ser atingido antes da precisão pedida. A coluna "Diferença" da tabela final, que compara com o método direto, mostra quando isso acontece.
- **Simpson** exige pontos igualmente espaçados e número par de subintervalos.
- **Interpolação** vai até o grau 4, como pede o enunciado.
- **Caracteres de moldura** (╔ ═ │): aparecem corretamente no Windows Terminal, no VS Code e em terminais Linux/macOS. Em consoles antigos do Windows, se aparecerem símbolos estranhos, use o Windows Terminal ou rode `chcp 65001` antes do programa.


Trabalho desenvolvido para a disciplina de Cálculo Numérico – UNIVASF, semestre 2026.2.
