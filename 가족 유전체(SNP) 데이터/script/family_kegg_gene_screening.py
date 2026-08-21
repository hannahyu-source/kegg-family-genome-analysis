"""가족 SNP × KEGG 질환연관 유전자(variant.txt 909개) 관심영역 스크리닝.

원래 아이디어: variant.txt 유전자를 GRCh37 좌표로 매핑해 가족 SNP 파일에서
해당 유전자 영역의 rsid를 직접 필터링하려 했다. 하지만 Phase 3에서 이미
NCBI ClinVar variant_summary.txt를 받아 rsid 단위로 매칭해뒀고, 그 파일에
GeneSymbol 컬럼이 그대로 들어있다 — 그래서 별도로 유전자 좌표를 받아올 필요 없이
family_clinvar_matches.csv를 그대로 재사용해서 "가족 SNP가 어느 유전자에 속하는지"를
알 수 있다. (ClinVar 자체가 GRCh37 좌표 기준으로 gene을 이미 매핑해뒀기 때문.)

이 스크립트가 하는 일:
  1. kegg_disease_gene_panel.csv(909개 유전자)와 family_clinvar_matches.csv를
     유전자 심볼로 join — "가족이 SNP를 가진 유전자 중 KEGG가 질환연관으로 등록한 것"
  2. 유전자별로 variant.txt 표현형 -> 연결 질병 -> 연결 network -> 연결 약물/약물군까지
     체인을 자동으로 이어붙인 리포트 생성 (network 이름은 kegg_variant.csv의
     network_ids를 kegg_network.csv와 join해서 보강)

출력 (가족 유전체(SNP) 데이터/output/):
  family_kegg_gene_screening.csv - 유전자별 스크리닝 결과 (가족 보유 여부 + 전체 체인)
"""

import sys
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
FAMILY_OUTPUT_DIR = SCRIPT_DIR.parent / "output"
KEGG_OUTPUT_DIR = SCRIPT_DIR.parent.parent / "KEGG 데이터" / "output"

MEMBER_ORDER = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]


def build_gene_network_names():
    """kegg_variant.csv의 network_ids를 kegg_network.csv 이름으로 풀어 유전자별로 모은다."""
    variant_df = pd.read_csv(KEGG_OUTPUT_DIR / "kegg_variant.csv", dtype=str).fillna("")
    network_df = pd.read_csv(KEGG_OUTPUT_DIR / "kegg_network.csv", dtype=str).fillna("")
    network_name = dict(zip(network_df["entry_id"], network_df["name"]))

    gene_networks = {}
    for gene, group in variant_df.groupby("gene_symbol"):
        ids = sorted({n for cell in group["network_ids"] for n in cell.split(";") if n})
        names = [network_name.get(n, n) for n in ids]
        gene_networks[gene] = "; ".join(names)
    return gene_networks


def main():
    panel = pd.read_csv(FAMILY_OUTPUT_DIR / "kegg_disease_gene_panel.csv", dtype=str).fillna("")
    matches = pd.read_csv(FAMILY_OUTPUT_DIR / "family_clinvar_matches.csv", dtype=str).fillna("")
    gene_networks = build_gene_network_names()

    panel_genes = set(panel["gene_symbol"])
    carrier_cols = [f"{m}_carries_alt_allele" for m in MEMBER_ORDER]

    # 이 유전자 영역에 "칩이 테스트한 rsid가 있다"와 "가족이 실제로 그 대립유전자를 가졌다"는
    # 전혀 다른 이야기다. BRCA1/BRCA2처럼 병원성 기록이 방대한 유전자는 테스트된 위치만으로도
    # 매칭이 많이 잡히지만, family_carries_anyone으로 걸러야 "실제 보유"만 남는다.
    matches["family_carries_anyone"] = (matches[carrier_cols] == "True").any(axis=1)

    family_hits_all = matches[matches["GeneSymbol"].isin(panel_genes)].copy()
    family_hits_carried = family_hits_all[family_hits_all["family_carries_anyone"]].copy()

    print(f"KEGG 질환연관 유전자(909개) 중 가족 SNP가 매칭된 유전자: {family_hits_all['GeneSymbol'].nunique()}개 (테스트된 위치 기준, 보유 여부 무관)")
    print(f"그중 가족이 실제로 대립유전자를 보유한 유전자: {family_hits_carried['GeneSymbol'].nunique()}개")
    print(f"실제 보유 rsid 총 {len(family_hits_carried)}건")

    panel_indexed = panel.set_index("gene_symbol")

    rows = []
    for gene, group in family_hits_carried.groupby("GeneSymbol"):
        panel_row = panel_indexed.loc[gene]

        carriers = set()
        for _, r in group.iterrows():
            for m, col in zip(MEMBER_ORDER, carrier_cols):
                if r[col] == "True":
                    carriers.add(m)

        significances = sorted(set(group["ClinicalSignificance"]) - {"", "-"})
        has_pathogenic = any("Pathogenic" in s and "Conflicting" not in s for s in significances)
        has_drug_response = any("drug response" in s for s in significances)

        rows.append(
            {
                "gene_symbol": gene,
                "kegg_variant_phenotypes": panel_row["variant_phenotypes"],
                "kegg_linked_diseases": panel_row["linked_disease_names"],
                "kegg_linked_networks": gene_networks.get(gene, ""),
                "kegg_linked_drug_targets": panel_row["linked_drug_target_names"],
                "family_carried_rsid_count": len(group),
                "family_carriers": ", ".join(sorted(carriers)),
                "clinvar_significances_carried": "; ".join(significances),
                "flag_pathogenic": has_pathogenic,
                "flag_drug_response": has_drug_response,
            }
        )

    result = pd.DataFrame(rows)
    result["priority"] = (
        result["flag_pathogenic"].astype(int) * 2
        + result["flag_drug_response"].astype(int)
        + (result["kegg_linked_diseases"].str.len() > 0).astype(int)
    )
    result = result.sort_values(["priority", "family_carried_rsid_count"], ascending=[False, False]).drop(columns="priority")
    result.to_csv(FAMILY_OUTPUT_DIR / "family_kegg_gene_screening.csv", index=False, encoding="utf-8-sig")

    print()
    print("=== 우선순위 상위 10개 유전자 ===")
    pd.set_option("display.max_colwidth", 35)
    print(
        result.head(10)[
            ["gene_symbol", "family_carriers", "flag_pathogenic", "flag_drug_response", "kegg_linked_diseases"]
        ].to_string(index=False)
    )
    print()
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_kegg_gene_screening.csv'}")


if __name__ == "__main__":
    main()
