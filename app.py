import re
from io import BytesIO
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import patheffects

# Opcional: Plotly para gráficos interativos
try:
    import plotly.express as px
    HAS_PLOTLY = True
except Exception:
    HAS_PLOTLY = False

st.set_page_config(page_title="Análise de Suscetibilidade", layout="wide")

# ---------------------------
# Constantes e utilitários
# ---------------------------
ABX_LEGENDA = [
    "GEN", "TOB", "AMI",
    "ATM",
    "CRO", "CAZ", "CTX", "CFO", "CPM",
    "AMC", "AMP", "PPT", "CZA",
    "MER", "IMP",
    "CIP", "LEV",
    "SUT",
    "POLI B",
]

def norm_str(s):
    if pd.isna(s):
        return ""
    return str(s).replace("–", "-").strip()

def extrai_num_codigo(cod):
    s = norm_str(cod).upper().replace(" ", "")
    m = re.search(r"\bMA\s*(\d+)", s)
    return int(m.group(1)) if m else np.nan

def normaliza_status(v):
    """Extrai o primeiro rótulo válido {SSD, R, S, I} em qualquer ponto da célula."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).upper().replace("–", "-")
    s = s.replace("INTERMEDIARIO", "I").replace("INTERMEDIÁRIO", "I")
    s = s.replace("SSD/I", "SSD").replace("SSD –", "SSD").replace("SSD–", "SSD")
    m = re.search(r"\b(SSD|R|S|I)\b", s)
    return m.group(1) if m else None

def achar_col(df, patt):
    for c in df.columns:
        if re.search(patt, norm_str(c), flags=re.IGNORECASE):
            return c
    return None

def perfil_sri_rows(df_in, abx_cols):
    """Calcula métricas por antibiótico.
    Retorna: lista de linhas [Antibiótico, N, Total, %Cobertura, R, S, I, SSD, %R, %S, %I+SSD, Sem Resultado]
    """
    rows = []
    total = len(df_in)
    for c in abx_cols:
        vals = [normaliza_status(v) for v in df_in[c]]
        valid = [v for v in vals if v is not None]
        n = len(valid)
        cnt = Counter(valid)
        S = cnt.get("S", 0)
        R = cnt.get("R", 0)
        I = cnt.get("I", 0)
        SSD = cnt.get("SSD", 0)
        perc_R = round(100 * R / n, 1) if n > 0 else 0.0
        perc_S = round(100 * S / n, 1) if n > 0 else 0.0
        perc_I_SSD = round(100 * (I + SSD) / n, 1) if n > 0 else 0.0
        cobertura = round(100 * n / total, 1) if total > 0 else 0.0
        rows.append([c, n, total, cobertura, R, S, I, SSD, perc_R, perc_S, perc_I_SSD, total - n])
    rows.sort(key=lambda r: (r[8], r[3], r[0]), reverse=True)
    return rows

def rows_to_df(rows):
    cols = [
        "Antibiótico", "N", "Total", "%Cobertura", "R", "S", "I", "SSD",
        "%R", "%S", "%I+SSD", "Sem Resultado",
    ]
    return pd.DataFrame(rows, columns=cols)

def format_df_numeric(df):
    fmt_cols = ["%Cobertura", "%R", "%S", "%I+SSD"]
    for c in fmt_cols:
        if c in df.columns:
            df[c] = df[c].astype(float).round(1)
    return df

def fig_topR(rows, titulo, top=15):
    dfp = rows_to_df(rows)
    if dfp.empty:
        return None
    dfp = dfp.sort_values("%R", ascending=False).head(top)
    if HAS_PLOTLY:
        fig = px.bar(
            dfp.sort_values("%R"), x="%R", y="Antibiótico", orientation="h",
            title=titulo, labels={"%R": "% Resistência (R)"}, text="%R"
        )
        # Mostrar rótulos de % dentro das barras
        fig.update_traces(texttemplate='%{x:.1f}% ', textposition='inside', textfont_color='white')
        fig.update_layout(yaxis={'categoryorder':'total ascending'})
        return fig
    else:
        fig, ax = plt.subplots(figsize=(6, max(3, len(dfp) * 0.35)))
        bars = ax.barh(dfp["Antibiótico"], dfp["%R"])
        ax.set_title(titulo)
        ax.set_xlabel("% Resistência (R)")
        # Rótulos dentro das barras
        try:
            ax.bar_label(bars, labels=[f"{v:.1f}%" for v in dfp["%R"]], label_type='center', color='white')
        except Exception:
            pass
        plt.tight_layout()
        return fig

def fig_stacked(rows, titulo, ordenar_por="%R", label_threshold=5.0, label_decimals=1, label_fontsize=11, label_vertical=True, label_vertical_dir="up", show_labels=True):
    dfp = rows_to_df(rows)
    if dfp.empty:
        return None
    if ordenar_por in dfp.columns:
        dfp = dfp.sort_values(ordenar_por, ascending=False)
    if HAS_PLOTLY:
        dfm = dfp.melt(id_vars=["Antibiótico", "N"], value_vars=["%S", "%I+SSD", "%R"],
                       var_name="Classe", value_name="Percentual")
        ordem = {"%S": "Sensível (S)", "%I+SSD": "Intermediário (I/SSD)", "%R": "Resistente (R)"}
        dfm["Classe"] = dfm["Classe"].map(ordem)
        # Monta rótulo condicional por ponto
        fmt = f"{{:.{label_decimals}f}}%"
        if show_labels and label_vertical:
            dfm["Label"] = dfm["Percentual"].apply(lambda v: fmt.format(float(v)) if float(v) >= float(label_threshold) else "")
        else:
            dfm["Label"] = ""
        fig = px.bar(
            dfm, x="Antibiótico", y="Percentual", color="Classe", title=titulo,
            color_discrete_map={
                "Sensível (S)": "#2ca02c",
                "Intermediário (I/SSD)": "#1f77b4",
                "Resistente (R)": "#d62728",
            }, text="Label"
        )
        # Rótulos dentro das barras empilhadas, verticais e com margem para evitar corte no topo
        # 90 = baixo->cima; -90 (ou 270) = cima->baixo
        angle = 0
        if show_labels and label_vertical:
            angle = 90 if str(label_vertical_dir).lower() in ("up", "cima", "baixo->cima", "bottom-up", "bottom_to_top") else -90
        fig.update_traces(
            textposition=('inside' if show_labels and label_vertical else 'none'),
            textangle=angle,
            textfont_color='white',
            textfont_size=label_fontsize,
            cliponaxis=False
        )
        fig.update_layout(
            xaxis_tickangle=-60,
            yaxis_range=[0, 102],
            uniformtext_minsize=10,
            uniformtext_mode='show',
            margin=dict(t=60)
        )
        return fig
    else:
        x = np.arange(len(dfp))
        fig, ax = plt.subplots(figsize=(max(8, len(dfp) * 0.8), 6.2))
        bars_S = ax.bar(x, dfp["%S"], label="Sensível (S)", color="#2ca02c")
        bars_I = ax.bar(x, dfp["%I+SSD"], bottom=dfp["%S"], label="Intermediário (I/SSD)", color="#1f77b4")
        bars_R = ax.bar(x, dfp["%R"], bottom=dfp["%S"] + dfp["%I+SSD"], label="Resistente (R)", color="#d62728")
        ax.set_title(titulo)
        ax.set_ylabel("% de isolados testados")
        ax.set_xticks(x)
        ax.set_xticklabels(dfp["Antibiótico"], rotation=60, ha="right")
        ax.set_ylim(0, 102)  # pequena folga para rótulos
        ax.legend(loc="upper right")
        # Rótulos centralizados em cada segmento
        if show_labels and label_vertical:
            try:
                fmt = f"{{:.{label_decimals}f}}%"
                def labels(vals):
                    return [fmt.format(float(v)) if float(v) >= float(label_threshold) else "" for v in vals]
                # 90 = baixo->cima; -90 = cima->baixo
                rot = 90 if str(label_vertical_dir).lower() in ("up", "cima", "baixo->cima", "bottom-up", "bottom_to_top") else -90
                texts = []
                texts += ax.bar_label(bars_S, labels=labels(dfp["%S"]), label_type='center', color='white', rotation=rot, padding=0, fontsize=label_fontsize)
                texts += ax.bar_label(bars_I, labels=labels(dfp["%I+SSD"]), label_type='center', color='white', rotation=rot, padding=0, fontsize=label_fontsize)
                texts += ax.bar_label(bars_R, labels=labels(dfp["%R"]), label_type='center', color='white', rotation=rot, padding=0, fontsize=label_fontsize)
                # contorno para legibilidade
                for t in texts:
                    if hasattr(t, 'set_path_effects'):
                        t.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])
            except Exception:
                pass
        plt.tight_layout()
        return fig

def to_csv_download(df):
    buffer = BytesIO()
    df.to_csv(buffer, index=False, encoding="utf-8-sig")
    buffer.seek(0)
    return buffer

# ---------------------------
# UI
# ---------------------------
st.title("Análise de Suscetibilidade de Isolados")
st.caption("Carregue o arquivo Excel e explore métricas por antibiótico com gráficos e tabelas interativas.")

with st.sidebar:
    st.header("Configurações")
    # Escolha da fonte de dados
    data_source = st.radio("Fonte dos dados", ["Arquivo padrão (Cepas.xlsx)", "Upload de arquivo"], index=0)
    up = None
    if data_source == "Upload de arquivo":
        up = st.file_uploader("Carregar Excel (.xlsx)", type=["xlsx"])
    header_row = st.number_input("Linha do cabeçalho (1 = primeira linha)", min_value=1, value=2)
    codigo_limite = st.number_input("Limite numérico p/ 2025 (ex.: MA >= 180)", min_value=0, value=180)
    top_n = st.slider("Top N por %R (gráfico horizontal)", 5, 30, 15)
    ordenar_por = st.selectbox("Ordenar gráfico empilhado por", ["%R", "%S", "%Cobertura", "Antibiótico"])
    min_n = st.slider("N mínimo testado (filtro)", 0, 100, 0)
    min_cob = st.slider("Cobertura mínima %", 0, 100, 0)
    usar_plotly = st.checkbox("Forçar gráficos interativos (Plotly)", value=HAS_PLOTLY)
    st.subheader("Rótulos das barras (empilhado)")
    show_labels = st.checkbox("Exibir rótulos internos (vertical)", value=True)
    label_threshold = st.slider("Mostrar rótulo a partir de (%)", 0, 20, 5, disabled=not show_labels)
    label_decimals = st.select_slider("Casas decimais", options=[0,1,2], value=1, disabled=not show_labels)
    label_fontsize = st.slider("Tamanho do texto", 8, 18, 11, disabled=not show_labels)
    label_vertical_dir = st.selectbox("Direção do rótulo vertical", ["Baixo -> cima", "Cima -> baixo"], index=0, disabled=not show_labels)
    st.divider()
    st.subheader("Espécies (MALDI-TOF)")
    especies_chart_tipo = st.selectbox("Tipo de gráfico", ["Barra", "Pizza"], index=0)
    especies_top_n = st.slider("Top N espécies", 3, 20, 10)
    especies_agrupar_outros = st.checkbox("Agrupar demais em 'Outros'", value=True)

# Define o arquivo a ser lido
DEFAULT_XLSX = Path("Cepas.xlsx")
file_to_read = None
source_label = ""

if data_source == "Arquivo padrão (Cepas.xlsx)":
    if DEFAULT_XLSX.exists():
        file_to_read = DEFAULT_XLSX
        source_label = f"Arquivo padrão: {DEFAULT_XLSX.name}"
    else:
        st.warning("Arquivo padrão 'Cepas.xlsx' não encontrado na pasta do app. Selecione 'Upload de arquivo'.")
        st.stop()
else:
    if up is None:
        st.info("Envie um arquivo Excel (.xlsx). O app espera cabeçalho na 2ª linha por padrão.")
        st.stop()
    file_to_read = up
    source_label = f"Upload: {getattr(up, 'name', 'arquivo')}"

# ---------------------------
# Leitura do Excel
# ---------------------------
try:
    df = pd.read_excel(file_to_read, dtype=object, header=int(header_row - 1))
except Exception as e:
    st.error(f"Erro ao ler a planilha: {e}")
    st.stop()

df.columns = [norm_str(c) for c in df.columns]
st.caption(f"Fonte dos dados: {source_label}")

col_codigo = achar_col(df, r"^CÓDIGO\s+UFPB$")
if not col_codigo:
    st.error("Coluna 'CÓDIGO UFPB' não encontrada. Verifique o cabeçalho correto.")
    st.dataframe(pd.DataFrame({"Colunas encontradas": df.columns}))
    st.stop()

col_maldi = achar_col(df, r"MALDLI-TOF|MALDI-TOF")
abx_cols_existentes = [c for c in ABX_LEGENDA if c in df.columns]
if not abx_cols_existentes:
    st.error("Nenhuma coluna de antibiótico da legenda foi encontrada no arquivo.")
    st.write("Colunas esperadas (exemplos):", ABX_LEGENDA)
    st.write("Colunas do arquivo:", list(df.columns))
    st.stop()

# Seleção dos antibióticos
st.sidebar.subheader("Antibióticos")
sel_abx = st.sidebar.multiselect("Selecionar colunas de antibióticos", options=abx_cols_existentes, default=abx_cols_existentes)
if not sel_abx:
    st.warning("Selecione pelo menos um antibiótico.")
    st.stop()

# Calcula partições
df["_NUM"] = df[col_codigo].apply(extrai_num_codigo)
grp_2025 = df[df["_NUM"] >= codigo_limite].copy()
grp_prev = df[(df["_NUM"] < codigo_limite) & (~df["_NUM"].isna())].copy()

tabs = st.tabs(["2025 / atual", "Anos anteriores", "Espécies (MALDI-TOF)"])

# Helper para uma aba de análise
def render_tab(grp_df, titulo):
    st.subheader(titulo)
    total_iso = len(grp_df)
    st.caption(f"Total de isolados no grupo: {total_iso}")
    if total_iso == 0:
        st.info("Sem registros neste grupo.")
        return

    rows = perfil_sri_rows(grp_df, sel_abx)
    df_rows = format_df_numeric(rows_to_df(rows))
    # Filtros
    df_rows = df_rows[(df_rows["N"] >= min_n) & (df_rows["%Cobertura"] >= min_cob)]
    if df_rows.empty:
        st.warning("Nenhum antibiótico atende aos filtros aplicados.")
        return

    c1, c2 = st.columns([2, 1])
    with c1:
        st.markdown("### Tabela de suscetibilidade")
        st.dataframe(
            df_rows.style
            .background_gradient(subset=["%R"], cmap="Reds")
            .background_gradient(subset=["%S"], cmap="Greens")
            .format({"%R": "{:.1f}", "%S": "{:.1f}", "%I+SSD": "{:.1f}", "%Cobertura": "{:.1f}"})
        , use_container_width=True)
    with c2:
        st.download_button(
            label="Baixar CSV",
            data=to_csv_download(df_rows),
            file_name=f"perfil_suscetibilidade_{titulo.replace(' ', '_').lower()}.csv",
            mime="text/csv",
        )

    c3, c4 = st.columns(2)
    with c3:
        st.markdown("### Top %R por antibiótico")
        fig1 = fig_topR(df_rows.values.tolist(), f"Top %R - {titulo}", top=top_n)
        if fig1 is not None:
            if HAS_PLOTLY and usar_plotly:
                st.plotly_chart(fig1, use_container_width=True)
            else:
                st.pyplot(fig1, use_container_width=True)
    with c4:
        st.markdown("### Suscetibilidade dos isolados (empilhado)")
        fig2 = fig_stacked(
            df_rows.values.tolist(), f"Suscetibilidade - {titulo}", ordenar_por=ordenar_por,
            label_threshold=label_threshold, label_decimals=label_decimals,
            label_fontsize=label_fontsize, label_vertical=show_labels,
            label_vertical_dir=("up" if label_vertical_dir.startswith("Baixo") else "down"),
            show_labels=show_labels
        )
        if fig2 is not None:
            if HAS_PLOTLY and usar_plotly:
                st.plotly_chart(fig2, use_container_width=True)
            else:
                st.pyplot(fig2, use_container_width=True)


with tabs[0]:
    render_tab(grp_2025, "2025 / atual")

with tabs[1]:
    render_tab(grp_prev, "Anos anteriores")

with tabs[2]:
    st.subheader("Distribuição de espécies (MALDI-TOF)")
    if col_maldi and col_maldi in df.columns:
        esp = df[col_maldi].dropna().astype(str).str.strip()
        esp = esp[(esp != "") & (esp != "*")]
        if not esp.empty:
            dist = esp.value_counts().sort_values(ascending=False)
            total = int(dist.sum())
            # Preparar Top N + 'Outros' se habilitado
            if especies_agrupar_outros and len(dist) > especies_top_n:
                top = dist.head(especies_top_n)
                outros_val = int(dist.iloc[especies_top_n:].sum())
                dist_plot = top.append(pd.Series({"Outros": outros_val}))
            else:
                dist_plot = dist
            dist_df = dist_plot.reset_index()
            dist_df.columns = ["Espécie", "Frequência"]
            dist_df["%"] = (dist_df["Frequência"] / dist_df["Frequência"].sum() * 100).round(1)

            st.dataframe(dist_df, use_container_width=True)

            if especies_chart_tipo in ("Pizza"):
                if HAS_PLOTLY and usar_plotly:
                    fig = px.pie(dist_df, values="Frequência", names="Espécie",
                                  hole=0.4, title="Distribuição de espécies (MALDI-TOF)")
                    fig.update_traces(textposition='inside', textinfo='percent+label')
                    # Aumenta levemente o tamanho do gráfico (altura); largura segue o container
                    fig.update_layout(height=600)
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    # Aumenta levemente o tamanho do gráfico em polegadas (largura x altura)
                    fig, ax = plt.subplots(figsize=(8, 6))
                    wedges, texts, autotexts = ax.pie(
                        dist_df["Frequência"], labels=dist_df["Espécie"], autopct='%1.1f%%',
                        startangle=90, wedgeprops=dict(width=0.4)
                    )
                    ax.set_title("Distribuição de espécies (MALDI-TOF)")
                    ax.axis('equal')
                    plt.tight_layout()
                    st.pyplot(fig, use_container_width=True)
            else:
                if HAS_PLOTLY and usar_plotly:
                    fig = px.bar(dist_df, x="Frequência", y="Espécie", orientation="h",
                                 title="Distribuição de espécies (MALDI-TOF)")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    fig, ax = plt.subplots()
                    ax.barh(dist_plot.index[::-1], dist_plot.values[::-1])
                    ax.set_title("Distribuição de espécies (MALDI-TOF)")
                    ax.set_xlabel("Frequência")
                    ax.set_ylabel("Espécie")
                    plt.tight_layout()
                    st.pyplot(fig, use_container_width=True)
        else:
            st.info("Sem valores válidos em MALDI-TOF.")
    else:
        st.warning("Coluna MALDI-TOF não encontrada; análise de espécies omitida.")

# ---------------------------
# Rodapé (guia único)
# ---------------------------
st.divider()
with st.expander("Guia rápido: uso, regras e cálculos", expanded=False):
        st.markdown(
                """
                • Como usar (barra lateral)
                - Passo 1 — Fonte dos dados
                    - Arquivo padrão (Cepas.xlsx): usa automaticamente o arquivo local se existir na pasta do app.
                    - Upload de .xlsx: envie outro arquivo para análise.
                - Passo 2 — Cabeçalho e leitura
                    - Linha do cabeçalho: 2 por padrão (porque a 1ª linha costuma ser título/legenda). Ajuste se a sua planilha tiver o cabeçalho em outra linha.
                - Passo 3 — Particionamento temporal
                    - Limite p/ 2025: separa os grupos a partir do número extraído de `CÓDIGO UFPB` após “MA”. Ex.: `MA180` ⇒ 180. Valores `>= limite` entram em “2025 / atual”; os demais, em “Anos anteriores”.
                - Passo 4 — Antibióticos e cobertura
                    - Antibióticos: selecione as colunas a considerar (somente as existentes na planilha aparecem).
                    - N mínimo testado: oculta antibióticos com menos de N resultados válidos.
                    - Cobertura mínima %: oculta antibióticos pouco testados no grupo (`N/Total`).
                - Passo 5 — Gráficos por antibiótico
                    - Top N por %R: define quantos antibióticos aparecem no gráfico horizontal de Top Resistência (%R).
                    - Ordenar gráfico empilhado por: ordena as barras por `%R`, `%S`, `%Cobertura` ou por nome.
                    - Rótulos das barras (empilhado): ative se quiser mostrar valores dentro das barras; ajuste limiar de exibição, casas decimais, tamanho e direção do texto.
                    - Forçar gráficos interativos (Plotly): liga/desliga versão interativa (útil para zoom e hover; desative se estiver pesado no navegador).
                - Passo 6 — Espécies (MALDI-TOF)
                    - Tipo de gráfico: “Barra” (horizontal) ou “Pizza”.
                    - Top N espécies: limita às mais frequentes.
                    - Agrupar demais em “Outros”: soma as demais espécies em um único item.
                - Dicas rápidas
                    - Se algum antibiótico “sumir”, verifique se os filtros de N mínimo/Cobertura mínima não estão altos demais.
                    - Se o Arquivo padrão não existir, use Upload.
                    - Use “Baixar CSV” para auditar e conferir casos específicos.

                • Como os dados são lidos
                - Normalização dos nomes de colunas (trim e troca de `–` por `-`).
                - Colunas‑chave: `CÓDIGO UFPB` (exato após normalização) e `MALDI-TOF` (aceita variações contendo “MALDI”).
                - Antibióticos considerados (usa apenas os que existem na planilha):
                    GEN, TOB, AMI, ATM, CRO, CAZ, CTX, CFO, CPM, AMC, AMP, PPT, CZA, MER, IMP, CIP, LEV, SUT, POLI B.

                • Regras e cálculos
                - `_NUM`: extraído de `CÓDIGO UFPB` com `MA(\d+)` (ex.: `MA24B`, `MA180`).
                - Particionamento: `_NUM >= 180` ⇒ “2025 / atual”; `_NUM < 180` ⇒ “Anos anteriores”.
                - Resultado por célula: pega o primeiro rótulo válido em `{SSD, R, S, I}`; aceita variações (`SSD/I`, `SSD–`, `R *`, `S (ok)`), e `INTERMEDIARIO/INTERMEDIÁRIO` vira `I`. Ignora vazio, `*`, `-`, `NaN`.
                - Para cada antibiótico: `N`, contagens `R/S/I/SSD`, `%R`, `%S`, `%I+SSD` (base `N`), `Total`, `%Cobertura = 100*N/Total`, `Sem Resultado = Total - N`.
                - Ordenação da tabela: `%R` (desc), depois `%Cobertura` (desc) e nome.

                • Saídas
                - Tabela: `Antibiótico | N | Total | %Cobertura | R | S | I | SSD | %R | %S | %I+SSD | Sem Resultado` (percentuais com 1 casa decimal) + botão “Baixar CSV”.
                - Gráficos: Top %R (horizontal) e gráfico empilhado de `%S`, `%I+SSD`, `%R` com rótulos opcionais.
                - Espécies (MALDI-TOF): remove `""` e `*`; mostra barras horizontais (ou pizza) com Top N e “Outros”.

                • Salvaguardas
                - Não eliminamos testes válidos: qualquer célula que contenha `SSD|R|S|I` é contada uma vez.
                - Diferenças de `N` refletem a cobertura real (células vazias/asteriscos não contam).
                - Conferência: exporte o CSV e verifique casos específicos (ex.: `GEN` em `>=180`).

                • Fórmulas
                - `%R = 100 * R / N`
                - `%S = 100 * S / N`
                - `%I+SSD = 100 * (I + SSD) / N`
                - `%Cobertura = 100 * N / Total`
                """
        )
