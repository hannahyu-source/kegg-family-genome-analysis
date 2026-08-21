"""results/figures/ 아래에 포트폴리오용 정적 그림(PNG) 4종을 생성한다.

파이프라인 개요 그림은 README.md의 Mermaid 다이어그램으로 대체하고(GitHub가 직접
렌더링하므로 별도 정적 이미지가 불필요), 여기서는 실제 분석 결과를 보여주는 그림만
만든다:
  kegg_relationship_summary.png - KEGG 관계 그래프 타입별 엣지 수 상위 10
  family_clinvar_summary.png    - 가족 ClinVar 매칭 결과의 임상적 유의성 분포(전 가족 합산)
  cfh_family_genotype.png       - CFH 변이 2건의 가족 구성원별 보유 여부
  pgx_family_heatmap.png        - 상위 PGx 유전자의 가족 구성원별 보유 변이 수

색상은 dataviz 스킬의 검증된 기본 팔레트(references/palette.md)를 그대로 사용한다:
크기(magnitude) 비교는 순차(sequential) 블루 1색, 상태(status) 그리드는 범주형 slot 1
(blue)만 "보유"를 나타내는 데 쓰고 나머지는 중립 회색 + 해칭으로 구분한다.
"""

import warnings
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results" / "tables"
PROCESSED_DIR = ROOT / "data" / "processed"
FIGURES_DIR = ROOT / "results" / "figures"

# dataviz 스킬 reference/palette.md 값 그대로
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"
BLUE = "#2a78d6"
BLUE_SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#0d366b"]
NEUTRAL_LIGHT = "#e3e2dd"
NEUTRAL_MID = "#9c9b95"

# 한글 라벨을 그리려면 한글 글리프가 있는 폰트가 필요하다(matplotlib 기본 폰트인
# DejaVu Sans는 한글을 지원하지 않아 네모(tofu)로 깨진다). OS별 대표 한글 폰트를
# 순서대로 시도하고, 하나도 없으면 경고만 남기고 기본 폰트로 진행한다.
_KOREAN_FONT_CANDIDATES = ["Malgun Gothic", "AppleGothic", "NanumGothic", "Noto Sans CJK KR"]
_available = {f.name for f in fm.fontManager.ttflist}
_korean_font = next((f for f in _KOREAN_FONT_CANDIDATES if f in _available), None)
if _korean_font is None:
    warnings.warn(
        "한글 지원 폰트를 찾지 못했습니다 (Malgun Gothic/AppleGothic/NanumGothic/Noto Sans CJK KR). "
        "그림의 한글 라벨이 깨질 수 있습니다. Linux에서는 `fonts-nanum` 패키지 설치를 권장합니다."
    )

plt.rcParams.update(
    {
        "font.family": _korean_font or "DejaVu Sans",
        "axes.unicode_minus": False,
        "font.size": 11,
        "text.color": TEXT_PRIMARY,
        "axes.edgecolor": GRID,
        "axes.labelcolor": TEXT_SECONDARY,
        "xtick.color": TEXT_SECONDARY,
        "ytick.color": TEXT_SECONDARY,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    }
)


def save(fig, name):
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"-> {path}")


def fig_kegg_relationship_summary():
    df = pd.read_csv(RESULTS_DIR / "kegg_graph_type_summary.csv")
    top = df.sort_values("edge_count", ascending=False).head(10).iloc[::-1]
    labels = [f"{s} → {t}" for s, t in zip(top["source_type"], top["target_type"])]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(labels, top["edge_count"], color=BLUE, height=0.6)
    for y, v in enumerate(top["edge_count"]):
        ax.text(v, y, f"  {v:,}", va="center", ha="left", color=TEXT_PRIMARY, fontsize=9)
    ax.set_xlabel("관계(edge) 수")
    ax.set_title("KEGG 엔트리 타입 간 관계 — 상위 10 (총 69,546건)", color=TEXT_PRIMARY, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    save(fig, "kegg_relationship_summary.png")


def fig_family_clinvar_summary():
    df = pd.read_csv(RESULTS_DIR / "family_clinvar_summary.csv")
    agg = df.groupby("clinical_significance")["variant_count"].sum().sort_values(ascending=False).head(10)
    agg = agg.iloc[::-1]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(agg.index, agg.values, color=BLUE, height=0.6)
    ax.set_xscale("log")
    for y, v in enumerate(agg.values):
        ax.text(v, y, f"  {v:,}", va="center", ha="left", color=TEXT_PRIMARY, fontsize=9)
    ax.set_xlabel("가족 5명 합산, alt allele 보유 매칭 건수 (로그 스케일)")
    ax.set_title("가족 ClinVar 매칭 결과 — 임상적 유의성 분포 상위 10", color=TEXT_PRIMARY, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, which="both", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    fig.text(
        0.01, -0.02,
        "참고: 대다수는 Benign(양성)이며, 병원성(Pathogenic) 계열은 소수다 — 매칭 = 병원성이 아님.",
        fontsize=8.5, color=TEXT_SECONDARY,
    )
    save(fig, "family_clinvar_summary.png")


def fig_cfh_family_genotype():
    members = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]
    variants = ["rs460897", "rs1061170"]
    matches = pd.read_csv(PROCESSED_DIR / "family_clinvar_matches.csv", dtype=str)
    cfh = matches[matches["rsid"].isin(variants)].drop_duplicates("rsid").set_index("rsid")

    fig, ax = plt.subplots(figsize=(7, 3))
    ax.set_xlim(0, len(members))
    ax.set_ylim(0, len(variants))
    ax.invert_yaxis()

    for row, variant in enumerate(variants):
        for col, member in enumerate(members):
            genotype = cfh.loc[variant, f"{member}_genotype"] if variant in cfh.index else None
            carries = cfh.loc[variant, f"{member}_carries_alt_allele"] if variant in cfh.index else None
            missing = not isinstance(genotype, str) or genotype in ("", "nan")

            if missing:
                face, edge, hatch = NEUTRAL_LIGHT, NEUTRAL_MID, "//"
                label = "측정 안 됨"
            elif carries == "True":
                face, edge, hatch = BLUE, BLUE, None
                label = genotype
            else:
                face, edge, hatch = "#ffffff", NEUTRAL_MID, None
                label = genotype

            ax.add_patch(
                Rectangle((col, row), 0.94, 0.86, facecolor=face, edgecolor=edge, hatch=hatch, linewidth=1.2)
            )
            text_color = "#ffffff" if (carries == "True" and not missing) else TEXT_PRIMARY
            ax.text(col + 0.47, row + 0.43, label, ha="center", va="center", fontsize=10, color=text_color)

    ax.set_xticks([c + 0.47 for c in range(len(members))])
    ax.set_xticklabels(members, color=TEXT_PRIMARY)
    ax.set_yticks([r + 0.43 for r in range(len(variants))])
    ax.set_yticklabels(variants, color=TEXT_PRIMARY)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("CFH 변이 2건 — 가족 구성원별 genotype", color=TEXT_PRIMARY, fontweight="bold", pad=14)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=BLUE, edgecolor=BLUE, label="대체 대립유전자 보유"),
        Rectangle((0, 0), 1, 1, facecolor="#ffffff", edgecolor=NEUTRAL_MID, label="참조형(비보유)"),
        Rectangle((0, 0), 1, 1, facecolor=NEUTRAL_LIGHT, edgecolor=NEUTRAL_MID, hatch="//", label="측정 안 됨"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False, fontsize=8.5)
    save(fig, "cfh_family_genotype.png")


def fig_pgx_family_heatmap():
    members = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]
    df = pd.read_csv(RESULTS_DIR / "family_pharmacogenomics.csv", dtype=str)
    top_genes = df["gene_symbol"].value_counts().head(12).index.tolist()

    counts = pd.DataFrame(0, index=top_genes, columns=members)
    for _, row in df[df["gene_symbol"].isin(top_genes)].iterrows():
        carriers = [c.strip() for c in row["carriers"].split(",") if c.strip()]
        for m in carriers:
            counts.loc[row["gene_symbol"], m] += 1

    vmax = max(counts.values.max(), 1)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list("blue_seq", BLUE_SEQ)
    im = ax.imshow(counts.values, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")

    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            v = counts.values[i, j]
            color = "#ffffff" if v / vmax > 0.55 else TEXT_PRIMARY
            ax.text(j, i, str(v), ha="center", va="center", fontsize=9, color=color)

    ax.set_xticks(range(len(members)))
    ax.set_xticklabels(members, rotation=30, ha="right", color=TEXT_PRIMARY)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes, color=TEXT_PRIMARY)
    ax.set_title("가족 PGx 변이 수 — 상위 12개 유전자 x 구성원", color=TEXT_PRIMARY, fontweight="bold", pad=12)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, shrink=0.7, label="보유 rsid 수 (drug response 계열)")
    cbar.outline.set_visible(False)
    save(fig, "pgx_family_heatmap.png")


def main():
    fig_kegg_relationship_summary()
    fig_family_clinvar_summary()
    fig_cfh_family_genotype()
    fig_pgx_family_heatmap()


if __name__ == "__main__":
    main()
