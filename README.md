# SMS com kernel de Bessel (primeiro tipo)

**GitHub:** [fabrizzioconde/SMS-Kernel-Bessel](https://github.com/fabrizzioconde/SMS-Kernel-Bessel)

Repositório de scripts em **R** para o método **SMS** (*Selection of Markers by Support vector machine*, seleção de marcadores com máquina de vetores de suporte em modo regressão — SVR), comparando kernels do SVR em dados **simulados** de genótipo e fenótipo (SNPs): **linear** e **radial (RBF / gaussiano)** via `e1071::svm`, e **várias configurações do kernel de Bessel** via `kernlab` (`besseldot`).

---

## Objetivo

Avaliar, em estudos de simulação com diferentes estruturas de efeitos e correlações entre marcadores, como o SMS recupera SNPs causais quando o SVR usa:

1. Kernel **linear** (`e1071::svm`);
2. Kernel **radial** — RBF / gaussiano (`e1071::svm`, `kernel = "radial"`, com parâmetro `gamma`), usualmente como **referência de kernel não linear clássico** frente ao Bessel;
3. Cinco variantes do kernel **Bessel** (`kernlab::ksvm` + `besseldot`), variando principalmente `sigma` e `order`.

O fluxo combina **importância da Random Forest**, **corte guiado pelo MSE do SVR** em grupos crescentes de SNPs e **refinamento por algoritmo genético (GA)**.

### Nota sobre o que está nos scripts `.R` desta pasta

A função `validacao_cruzada` está documentada no código como implementação para **kernel radial e linear** (`e1071`). Nos arquivos `r/grupo_*/SMS_Completo_Corr_*_grupo_*/SMS_Completo_Corr_*_Bessel_grupo_*.R`, o primeiro bloco do SMS (`i = 1`) está parametrizado com **`kernel = "linear"`** (corte, GA e saídas como `SMS Linear`). Para **reproduzir ou manter** também o ramo **radial** na mesma pipeline (corte + GA + `sink`), basta repetir esse bloco trocando para `kernel = "radial"` e definindo `gamma` (no código já existem valores de `gamma` em outros trechos, p.ex. na validação só com SNPs causais, que podem servir de ponto de partida).

---

## Estrutura do repositório

| Caminho | Conteúdo |
|---------|----------|
| `docs/thesis/` | **TCC** em PDF (`TCC_Diogo_Moura_Ferreira.pdf`) — trabalho de conclusão associado ao projeto. |
| `r/grupo_1/`, `r/grupo_2/`, `r/grupo_3/` | Três **réplicas** experimentais; em cada uma, subpastas `SMS_Completo_Corr_N_grupo_G/` por cenário. |
| `SMS_Completo_Corr_N_grupo_G/` (dentro de `r/grupo_G/`) | Um **cenário de simulação** `N` (1–6): scripts `.R`, saídas `.txt`, PDFs de diagnóstico e objetos `.RData`. |
| `pyproject.toml`, `uv.lock` | Projeto **Python** gerenciado por [uv](https://github.com/astral-sh/uv) (ambiente virtual e dependências; ver secção abaixo). |

Em cada pasta de cenário você encontra, em geral:

- **`SMS_Completo_Corr_N_Bessel_grupo_G.R`** — script principal: simulação, SMS completo, saídas e gráficos (em `r/grupo_G/SMS_Completo_Corr_N_grupo_G/`).
- **`SNPs_corte_selecionados_Corr_N_grupo_G.R`** (quando existir) — utilitário que cruza SNPs selecionados na **etapa de corte** (`snps_selec_corte`) com a lista de SNPs **causais** definida no cenário.
- **`SimulacaoN.RData`**, **`SMS_corr_N.RData`** — objetos R salvos após execução.
- **PDFs**: histograma e boxplot do fenótipo; MSE do SVR vs. número de marcadores; evolução do GA (ggplot2).
- **`my_output_SMS_Bessel_Corr_N.txt`** — texto com SNPs finais por kernel (nesta árvore: **SMS Linear** + Bessel 1–5; se você adicionar o ramo radial, incluir algo como `SMS Radial`), união, interseção, comparações com valor-*p* e métricas no conjunto **causal** conhecido.
- **`my_output_SMS_Rank_Corr_N.txt`** — saídas relacionadas a seleção por valor-*p* e rank da RF.
- **`rank_Random_Forest_N.txt`**, **`rank_global_N.txt`** — exportação do rank de importância da RF e tabelas com valor-*p* (formato CSV via `write.csv`).
- **`.Rhistory`** — histórico de comandos do R (gerado pelo ambiente local).

---

## Cenários de simulação (`Corr_1` … `Corr_6`)

Todos usam `scrime::simulateSNPglm` para gerar genótipo e fenótipo. Os **cenários diferem** em tamanho amostral, número de SNPs, posição dos causais, estrutura de interações (`list.ia`), `beta0` e `beta`.

| Cenário | Indivíduos | SNPs | Ideia resumida |
|---------|------------|------|----------------|
| **Corr_1** | 1000 | 100 | Oito SNPs causais (1–8), várias interações; `beta0=640`, efeitos unitários. |
| **Corr_2** | 1000 | 100 | Quatro grupos de dois SNPs; `beta=c(2,2,2,2)`. |
| **Corr_3** | 1000 | 100 | Três grupos de três SNPs; `beta=c(3,3,3)`. |
| **Corr_4** | 1000 | 100 | Estrutura mista (SNPs isolados e grupos); `beta` heterogêneo. |
| **Corr_5** | 1000 | 100 | Um bloco causal com quatro SNPs; `beta=c(4)`. |
| **Corr_6** | 250 | 1000 | Painel maior, sete causais nas posições 1, 10, …, 60; efeitos majoritariamente aditivos; `beta0=0`, um efeito maior (`900` no 4.º causal). |

Detalhes exatos estão nos comentários e parâmetros no início de cada `SMS_Completo_Corr_N_Bessel_grupo_G.R`.

---

## Ambiente virtual Python (uv)

O núcleo da análise é **R**; o **Python** com [uv](https://docs.astral.sh/uv/) fornece um console R moderno (**radian**) e **IPython** para tarefas auxiliares.

**Requisitos:** Python 3.11+ e [uv](https://docs.astral.sh/uv/getting-started/installation/) instalado (ou `python -m pip install uv`).

```bash
cd /caminho/para/SMS-Kernel-Bessel
uv sync --all-groups    # cria .venv e instala radian + ipython (grupo dev)
```

Ativar o ambiente e abrir o R via radian (com R instalado no sistema):

```bash
# Windows (cmd): .venv\Scripts\activate
# Git Bash: source .venv/Scripts/activate
radian
```

Dependências principais: `radian` (runtime); `ipython` (grupo `dev`, opcional).

---

## Fluxo do método (script principal)

1. **Pacotes**: `scrime`, `e1071`, `kernlab`, `randomForest`, `doParallel`, `GA`, `ggplot2` (instalação automática se faltar).
2. **Diretório de trabalho**: `setwd(dirname(rstudioapi::getActiveDocumentContext()$path))` — pressupõe execução no **RStudio** com o script como documento ativo.
3. **Simulação**: gera `dados[[1]]` (genótipo + `fenotipo`), teste de Shapiro-Wilk, histogramas e boxplots.
4. **Funções auxiliares**:
   - `validacao_cruzada` — SVR com `e1071::svm`; o **mesmo** laço de validação cruzada serve para **`kernel = "linear"`** ou **`kernel = "radial"`** (RBF), passando `gamma`, `cost` e `epsilon` como hoje (vide comentário “Kernel Radial e Linear” no `.R`).
   - `validacao_cruzada_bessel` — SVR com `kernlab::ksvm` e kernel `besseldot(sigma, order, degree)`.
   - `valor.p` — regressão linear SNP a SNP; valor-*p* bruto e ajuste tipo **Bonferroni** (`m * p`).
5. **Random Forest** (`randomForest`): ordenação dos SNPs por importância (`%IncMSE`).
6. **Corte (SMS)**: para cada kernel, acrescenta SNPs de 10 em 10 segundo o rank da RF; em cada passo calcula MSE médio da validação cruzada (10 folds). O **corte** usa o índice do **mínimo MSE** na curva (comentários no código mencionam comportamento no primeiro marcador).
7. **GA** (`GA::ga`, codificação binária): segunda etapa sobre o subconjunto pós-corte; **função de aptidão** = média da **correlação de Pearson** (predito vs. fenótipo) na validação cruzada, com o mesmo SVR/kernel da etapa correspondente.
8. **Agregação**: `uniao_snps` e `intersecao_snps` sobre as listas de SNPs escolhidos pelo GA para cada kernel.
9. **Referência “causal”**: no final, o script fixa nomes como `SNP1` … `SNP8` (ou o conjunto adequado ao cenário) e reporta correlação, R² ajustado e MSE do SVR **só com os causais verdadeiros**, para comparação com o SMS.

---

## Kernel radial (RBF) em relação ao Bessel

O **radial** em `e1071` é o kernel gaussiano \(K(x,x') = \exp(-\gamma\|x-x'\|^2)\). No SMS, ele usa a **mesma** estrutura que o linear: rank da RF → grupos de SNPs → mínimo do MSE na CV → GA maximizando correlação média na CV. Assim, a comparação **radial vs. Bessel** isola o efeito da *família* de kernel (RBF no `e1071` *versus* Bessel no `kernlab`), mantendo o restante do método idêntico.

---

## Kernels Bessel utilizados

No código, após o bloco com `e1071` no SMS (`i=1`, hoje **linear**), os índices `i=2…6` correspondem a **Bessel 1–5**, com parâmetros distintos, por exemplo:

- Bessel 1: `sigma=0.10`, `order=1`, `degree=1`
- Bessel 2: `sigma=0.30`, `order=1`, `degree=1`
- Bessel 3: `sigma=0.50`, `order=1`, `degree=1`
- Bessel 4: `sigma=0.50`, `order=1.20`, `degree=1`
- Bessel 5: `sigma=0.10`, `order=1.50`, `degree=1`

(`cost`, `epsilon` e `folds` seguem valores definidos em cada bloco.)

---

## Script `SNPs_corte_selecionados_*.R`

Define `selecionar_snp_causais(selected_snps, causal_snps)` e, em loop para `i` em 1…6, compara `snps_selec_corte[[i]]` com um vetor fixo de SNPs causais do cenário (ex.: `SNP1`, `SNP10`, … em Corr_6). **Deve ser executado no mesmo ambiente** em que o script principal já rodou (objetos `snps_selec_corte` na sessão).

---

## Como executar

1. Instalar [R](https://www.r-project.org/) e, de preferência, [RStudio](https://posit.co/download/rstudio-desktop/).
2. Abrir o arquivo `r/grupo_G/SMS_Completo_Corr_N_grupo_G/SMS_Completo_Corr_N_Bessel_grupo_G.R` correspondente ao cenário e ao grupo desejados.
3. Executar o script completo (pode levar bastante tempo por causa da RF, dos loops de SVR e do GA com `parallel=TRUE`).
4. Conferir PDFs e arquivos `.txt` na **mesma pasta** do script.

Para uso fora do RStudio, substitua o `setwd` por um caminho fixo para a pasta do script ou use `here::here()` / argumentos de linha de comando, pois `rstudioapi::getActiveDocumentContext()` exige o RStudio.

---

## Parâmetros globais recorrentes (exemplos)

- Random Forest: `ntree=4000`, `mtry = ncol(dados)-1`, `importance=TRUE`.
- Corte sequencial: `passo=10`, `percentual_snps=0.95` (limite superior de SNPs considerados no gráfico de MSE).
- GA (típico): `popSize=100`, `maxiter=10`, `run=30`, `pcross=0.8`, `pmut=0.1`, `elitism=5`.

---

## Resumo

Este projeto documenta uma **bateria de simulações** para o **SMS com SVR**, com ênfase no **kernel de Bessel (primeiro tipo)** em comparação com baselines **`e1071` (linear e, na mesma função, radial/RBF)** e com **três grupos de replicação** e **seis configurações de cenário** (`Corr_1`–`Corr_6`), mais análises auxiliares de SNPs no corte e saídas reprodutíveis em texto e gráficos.
