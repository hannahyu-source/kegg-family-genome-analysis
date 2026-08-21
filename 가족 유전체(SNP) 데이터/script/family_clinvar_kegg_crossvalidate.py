"""ClinVar 매칭 결과와 KEGG 질병연관 유전자 패널을 유전자 심볼로 join해서
"ClinVar가 병원성이라 하고 KEGG도 독립적으로 질병 연관을 인정하는" 항목만 걸러낸다.

두 데이터베이스는 서로 다른 큐레이션 과정을 거친 독립적인 소스다. 같은 유전자를
양쪽 다 질병과 연결짓고 있다면, ClinVar 단독 저신뢰(제출자 1명·무검토) 항목보다는
신뢰도가 한 단계 높다고 볼 수 있다 — 다만 "같은 유전자"이지 "같은 변이/기전"이라는
뜻은 아니라는 점은 여전히 유의해야 한다(아래 출력에도 명시).

필터 조건:
  1. family_clinvar_matches.csv에서 ClinicalSignificance가 Pathogenic 계열이고
     (Conflicting 제외) 실제로 가족 중 누군가 그 대립유전자를 보유
  2. 그 유전자가 kegg_disease_gene_panel.csv에 있고, linked_disease_count > 0
     (KEGG variant.txt 엔트리가 실제 DISEASE: H##### 필드로 질병을 명시한 경우만 인정 —
     단순히 variant.txt에 이름만 올라 있는 것은 불충분한 근거로 보고 제외)

출력 (가족 유전체(SNP) 데이터/output/):
  family_clinvar_kegg_crossvalidated.csv
"""

import sys
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
FAMILY_OUTPUT_DIR = SCRIPT_DIR.parent / "output"

MEMBER_ORDER = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]


def main():
    matches = pd.read_csv(FAMILY_OUTPUT_DIR / "family_clinvar_matches.csv", dtype=str)
    panel = pd.read_csv(FAMILY_OUTPUT_DIR / "kegg_disease_gene_panel.csv", dtype=str)

    # 1) ClinVar 쪽: Pathogenic 계열(Conflicting 제외) + 실제 보유자 존재
    is_pathogenic = matches["ClinicalSignificance"].str.contains("Pathogenic", na=False) & \
        ~matches["ClinicalSignificance"].str.contains("Conflicting", na=False)

    carrier_cols = [f"{m}_carries_alt_allele" for m in MEMBER_ORDER]
    carried_by_anyone = (matches[carrier_cols] == "True").any(axis=1)

    candidates = matches[is_pathogenic & carried_by_anyone].copy()
    print(f"ClinVar 병원성(비상충) + 실제 보유: {len(candidates)}건, 유전자 {candidates['GeneSymbol'].nunique()}개")

    # 2) KEGG 쪽: variant.txt에 명시적 DISEASE: H##### 링크가 있는 유전자만 인정
    panel["linked_disease_count"] = pd.to_numeric(panel["linked_disease_count"], errors="coerce").fillna(0)
    kegg_confirmed_genes = panel[panel["linked_disease_count"] > 0].set_index("gene_symbol")
    print(f"KEGG variant.txt에서 실제 질병(H#####)과 연결된 유전자: {len(kegg_confirmed_genes)}개 (전체 패널 {len(panel)}개 중)")

    # 3) join
    crossvalidated = candidates[candidates["GeneSymbol"].isin(kegg_confirmed_genes.index)].copy()
    crossvalidated["kegg_linked_disease_count"] = crossvalidated["GeneSymbol"].map(kegg_confirmed_genes["linked_disease_count"])
    crossvalidated["kegg_linked_disease_names"] = crossvalidated["GeneSymbol"].map(kegg_confirmed_genes["linked_disease_names"])
    crossvalidated["kegg_variant_phenotypes"] = crossvalidated["GeneSymbol"].map(kegg_confirmed_genes["variant_phenotypes"])

    carriers_summary = []
    for _, row in crossvalidated.iterrows():
        carriers = [m for m in MEMBER_ORDER if row[f"{m}_carries_alt_allele"] == "True"]
        carriers_summary.append(", ".join(carriers))
    crossvalidated["carriers"] = carriers_summary

    out_cols = [
        "rsid", "GeneSymbol", "ClinicalSignificance", "ReviewStatus", "NumberSubmitters",
        "PhenotypeList", "carriers",
        "kegg_linked_disease_count", "kegg_linked_disease_names", "kegg_variant_phenotypes",
    ] + [f"{m}_genotype" for m in MEMBER_ORDER]

    crossvalidated = crossvalidated[out_cols].sort_values("NumberSubmitters", ascending=False)
    crossvalidated.to_csv(FAMILY_OUTPUT_DIR / "family_clinvar_kegg_crossvalidated.csv", index=False, encoding="utf-8-sig")

    print()
    print(f"교차검증 통과: {len(crossvalidated)}건 (유전자 {crossvalidated['GeneSymbol'].nunique()}개)")
    if len(crossvalidated):
        pd.set_option("display.max_colwidth", 45)
        print(crossvalidated[["rsid", "GeneSymbol", "ClinicalSignificance", "NumberSubmitters", "carriers", "kegg_linked_disease_names"]].to_string(index=False))
    print()
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_clinvar_kegg_crossvalidated.csv'}")
    print()
    print("주의: '같은 유전자'가 두 DB에서 질병과 연결된다는 뜻이지, ClinVar가 병원성으로 본")
    print("바로 그 변이 기전을 KEGG가 확인해준다는 뜻은 아닙니다. 여전히 참고용입니다.")


if __name__ == "__main__":
    main()
