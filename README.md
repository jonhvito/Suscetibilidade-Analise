# Análise de Suscetibilidade — Documentação técnica

Aplicativo Streamlit para análise de resultados de suscetibilidade: leitura de Excel, padronização de resultados (R/S/I/SSD), métricas por antibiótico, gráficos (Top %R e empilhado) e distribuição de espécies (MALDI‑TOF), com exportação em CSV.

## Stack e dependências

- Python 3.8+
- Bibliotecas (ver versões em `requirements.txt`):
  - pandas, numpy, openpyxl
  - streamlit
  - matplotlib
  - plotly (opcional, para gráficos interativos)

Instalação:

```
pip install -r requirements.txt
```

Execução:

```
streamlit run app.py
```

## Formato de entrada

- Arquivo Excel: `Cepas.xlsx` (ou Upload via UI).
- Cabeçalho na 2ª linha (1ª linha é título/legenda). No aplicativo, isso corresponde a `header=1` (ajustável na barra lateral).
- Normalização de nomes de colunas: trim e substituição de `–` (en-dash) por `-` (hífen).
- Colunas‑chave esperadas:
  - `CÓDIGO UFPB` (detecção exata após normalização)
  - `MALDI-TOF` (ou variações contendo “MALDI”, p.ex. `MALDLI-TOF`)
- Antibióticos suportados (usa somente os que existirem na planilha):
  `GEN, TOB, AMI, ATM, CRO, CAZ, CTX, CFO, CPM, AMC, AMP, PPT, CZA, MER, IMP, CIP, LEV, SUT, POLI B`.

## Pipeline de processamento

1) Entrada e leitura
	- Carrega Excel com `pandas.read_excel(..., header=1)`. Os nomes das colunas são normalizados (`–` → `-`, trim).

2) Colunas‑chave
	- Identificador: `CÓDIGO UFPB`.
	- Espécies: coluna que contenha “MALDI” (`MALDLI-TOF` no arquivo original é suportada).
	- Antibióticos: lista fixa acima; somente colunas existentes são consideradas.

3) Extração do número da amostra
	- A partir de `CÓDIGO UFPB`, regex `MA(\d+)` (case‑insensitive) captura o inteiro em códigos como `MA24B`, `MA180`.
	- Campo derivado `_NUM` recebe o número extraído (ou `NaN` se não encontrado).

4) Particionamento temporal
	- Grupo “2025 / atual”: `_NUM >= 180` (inclusivo; `MA180` entra aqui).
	- Grupo “Anos anteriores”: `_NUM < 180` (somente quando `_NUM` não é nulo).

5) Padronização de resultados por antibiótico
	- Normalização robusta que extrai o primeiro rótulo válido por regex dentre `{SSD, R, S, I}`.
	- Tolera variações como `SSD/I`, `SSD –`, `SSD–`, `R *`, `S (ok)` e converte `INTERMEDIARIO/INTERMEDIÁRIO` em `I`.
	- Valores inválidos ignorados: vazio, `*`, `-`, `NaN`.

6) Métricas por antibiótico e grupo
	- Para cada coluna de antibiótico `c` dentro do grupo:
	  - `N`: nº de amostras com resultado válido em `c`.
	  - Contagens: `R`, `S`, `I`, `SSD`.
	  - Percentuais (base `N`): `%R`, `%S`, `%I+SSD`.
	  - Cobertura: `Total` = nº de amostras no grupo; `%Cobertura = 100*N/Total`; `Sem Resultado = Total - N`.
	- Ordenação: por `%R` (desc), depois por `%Cobertura` (desc) e, por fim, por nome.

7) Saída tabular
	- Cabeçalhos: `Antibiótico | N | Total | %Cobertura | R | S | I | SSD | %R | %S | %I+SSD | Sem Resultado`.
	- Percentuais formatados com 1 casa decimal.
	- Exportação CSV disponível via botão na interface.

8) Gráficos
	- Top %R: gráfico de barras horizontal (Top N por `%R`).
	- Empilhado: barras empilhadas de `%S`, `%I+SSD`, `%R`, com rótulos internos opcionais, orientação vertical configurável e limiar mínimo para exibição.

9) Distribuição de espécies (MALDI‑TOF)
	- Limpeza de valores `""` e `*`.
	- Barra horizontal ou pizza (opcional), com Top N e agrupamento de demais em “Outros”.

10) Salvaguardas
	- Não há filtro que elimine testes válidos: qualquer célula contendo `SSD|R|S|I` (mesmo com texto extra) conta uma única vez.
	- Diferenças de `N` por antibiótico refletem a cobertura real da planilha (células vazias/asteriscos não contam).
	- Conferência pontual: exporte a tabela para CSV e verifique casos específicos (ex.: `GEN` no grupo `_NUM >= 180`) cruzando `CÓDIGO UFPB` com o valor bruto.

### Fórmulas

- `%R = 100 * R / N`
- `%S = 100 * S / N`
- `%I+SSD = 100 * (I + SSD) / N`
- `%Cobertura = 100 * N / Total`

## Controles (barra lateral)

- Fonte dos dados: Arquivo padrão (Cepas.xlsx) ou Upload de .xlsx.
- Linha do cabeçalho: padrão 2.
- Limite p/ 2025: padrão 180; define o corte (`_NUM >= limite`).
- Antibióticos: multiseleção dentre as colunas existentes.
- N mínimo testado: padrão 0; filtra antibióticos com poucos resultados válidos.
- Cobertura mínima %: padrão 0; filtra antibióticos pouco testados.
- Top N por %R: padrão 15; controla o Top do gráfico horizontal.
- Ordenar gráfico empilhado por: `%R`, `%S`, `%Cobertura` ou `Antibiótico`.
- Rótulos (empilhado): habilitar/ocultar; limiar (%), casas decimais, tamanho, direção (baixo→cima/cima→baixo).
- Forçar gráficos interativos (Plotly): liga/desliga versão interativa.
- Espécies: tipo (Barra/Pizza), Top N (padrão 10) e “Agrupar demais em 'Outros'” (padrão ativado).

## Fluxo de uso sugerido

1) Selecione a fonte (Arquivo padrão ou Upload) e confirme a linha do cabeçalho.
2) Ajuste o limite p/ 2025 se necessário (default 180).
3) Escolha os antibióticos a analisar e aplique filtros de N mínimo/Cobertura mínima.
4) Examine a tabela; baixe o CSV para auditoria se precisar.
5) Explore os gráficos (Top %R e Empilhado) e a aba de espécies.

## Solução de problemas

- “Coluna 'CÓDIGO UFPB' não encontrada”: verifique a linha do cabeçalho e o nome exato após normalização (hifens `–` vs `-`). A UI exibe as colunas lidas para diagnóstico.
- “Nenhuma coluna de antibiótico encontrada”: confirme se a planilha usa os nomes esperados e se estão na segunda linha (cabeçalho).
- “Sem registros no grupo”: pode indicar que `_NUM` não foi extraído; confira o padrão `MA<numero>` na coluna `CÓDIGO UFPB`.
- Desempenho/Renderização: desative Plotly caso esteja pesado; reduza Top N ou o número de antibióticos selecionados.

## Observações

- O app não altera a planilha de origem; todas as transformações ocorrem em memória.
- O CSV baixado reflete exatamente os filtros e ordenações aplicados na interface.
