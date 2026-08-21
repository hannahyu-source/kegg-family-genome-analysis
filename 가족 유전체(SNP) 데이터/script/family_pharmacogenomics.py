"""가족 SNP 중 ClinVar가 "drug response"(약물반응)로 분류한 변이를 뽑아,
KEGG 약물군(dgroup)·약물(drug) 데이터와 연결한다.

이 방향은 Phase 3(가족 ClinVar 매칭)에서 의도적으로 손대지 않았던 categoy다 —
"Pathogenic"은 질병 위험이라 신중해야 했지만, "drug response"는 질병 진단이 아니라
약물 대사 속도/용량 반응에 관한 것이라 훨씬 낮은 민감도로 다룰 수 있다. 그래도 이 역시
참고용이며 실제 처방 변경은 의사·약사와 상의해야 한다는 점은 동일하다.

두 가지 방식으로 KEGG와 연결한다:
  1. 유전자 심볼이 KEGG dgroup 이름에 등장하는 경우 (예: "CYP2D6 substrate") —
     그 dgroup에 속한 모든 약물이 넓게 영향받을 수 있는 후보군
  2. ClinVar PhenotypeList에 적힌 "약물명 response" 패턴에서 약물명을 뽑아
     KEGG drug 엔트리와 직접 매칭 — 훨씬 좁고 구체적인 연결

출력 (가족 유전체(SNP) 데이터/output/):
  family_pharmacogenomics.csv         - 변이별 상세 (유전자, 보유자, 매칭된 KEGG dgroup/drug)
  family_pharmacogenomics_summary.csv - 구성원별 PGx 변이 수 · 영향권 KEGG 약물 수
"""

import re
import sys
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
FAMILY_OUTPUT_DIR = SCRIPT_DIR.parent / "output"
KEGG_OUTPUT_DIR = SCRIPT_DIR.parent.parent / "KEGG 데이터" / "output"

MEMBER_ORDER = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]
DRUG_NAME_PATTERN = re.compile(r"([A-Za-z][A-Za-z0-9/\-]{2,})\s+response\b")


def load_kegg_tables():
    dgroup_df = pd.read_csv(KEGG_OUTPUT_DIR / "kegg_dgroup.csv", dtype=str).fillna("")
    drug_df = pd.read_csv(KEGG_OUTPUT_DIR / "kegg_drug.csv", dtype=str).fillna("")
    drug_df["primary_name"] = drug_df["name"].str.split(";").str[0].str.split("(").str[0].str.strip()
    return dgroup_df, drug_df


def find_dgroups_for_gene(gene_symbol: str, dgroup_df: pd.DataFrame):
    mask = dgroup_df["name"].str.contains(rf"\b{re.escape(gene_symbol)}\b", case=False, na=False, regex=True)
    hits = dgroup_df[mask]
    return [
        f"{entry_id}:{name}({member_count}개 약물)"
        for entry_id, name, member_count in zip(hits["entry_id"], hits["name"], hits["member_count"])
    ]


def find_named_drugs(phenotype_list: str, drug_df: pd.DataFrame):
    candidates = {m.group(1).lower() for m in DRUG_NAME_PATTERN.finditer(phenotype_list)}
    if not candidates:
        return []

    name_lookup = dict(zip(drug_df["primary_name"].str.lower(), zip(drug_df["entry_id"], drug_df["primary_name"])))
    matched = []
    for cand in candidates:
        if cand in name_lookup:
            entry_id, primary_name = name_lookup[cand]
            matched.append(f"{entry_id}:{primary_name}")
    return sorted(set(matched))


def main():
    matches = pd.read_csv(FAMILY_OUTPUT_DIR / "family_clinvar_matches.csv", dtype=str).fillna("")
    dgroup_df, drug_df = load_kegg_tables()

    carrier_cols = [f"{m}_carries_alt_allele" for m in MEMBER_ORDER]
    carried_by_anyone = (matches[carrier_cols] == "True").any(axis=1)
    is_drug_response = matches["ClinicalSignificance"].str.contains("drug response", na=False)

    pgx = matches[is_drug_response & carried_by_anyone].copy()
    print(f"drug response 계열 변이 중 가족이 실제 보유: {len(pgx)}건, 유전자 {pgx['GeneSymbol'].nunique()}개")

    rows = []
    gene_dgroup_cache = {}
    for _, row in pgx.iterrows():
        gene = row["GeneSymbol"]
        if gene not in gene_dgroup_cache:
            gene_dgroup_cache[gene] = find_dgroups_for_gene(gene, dgroup_df)
        dgroups = gene_dgroup_cache[gene]

        named_drugs = find_named_drugs(row["PhenotypeList"], drug_df)
        carriers = [m for m in MEMBER_ORDER if row[f"{m}_carries_alt_allele"] == "True"]

        rows.append(
            {
                "rsid": row["rsid"],
                "gene_symbol": gene,
                "clinical_significance": row["ClinicalSignificance"],
                "phenotype_list": row["PhenotypeList"],
                "carriers": ", ".join(carriers),
                "matched_kegg_dgroups": "; ".join(dgroups),
                "matched_kegg_drugs": "; ".join(named_drugs),
            }
        )

    out_df = pd.DataFrame(rows).sort_values(["gene_symbol", "rsid"])
    out_df.to_csv(FAMILY_OUTPUT_DIR / "family_pharmacogenomics.csv", index=False, encoding="utf-8-sig")

    # 구성원별 요약: PGx 변이 수, 영향권 dgroup 수(중복제거), 그 dgroup들이 커버하는 약물 총합(중복제거)
    dgroup_member_ids = dict(zip(dgroup_df["entry_id"], dgroup_df["member_ids"]))
    summary_rows = []
    for member in MEMBER_ORDER:
        carrier_col = f"{member}_carries_alt_allele"
        member_variants = pgx[pgx[carrier_col] == "True"]
        genes = set(member_variants["GeneSymbol"])

        affected_dgroup_ids = set()
        for gene in genes:
            for entry in gene_dgroup_cache.get(gene, []):
                affected_dgroup_ids.add(entry.split(":")[0])

        affected_drug_ids = set()
        for did in affected_dgroup_ids:
            member_ids = dgroup_member_ids.get(did, "")
            affected_drug_ids.update(x for x in member_ids.split(";") if x)

        summary_rows.append(
            {
                "member": member,
                "pgx_variant_count": len(member_variants),
                "pgx_gene_count": len(genes),
                "pgx_genes": ", ".join(sorted(genes)),
                "affected_kegg_dgroup_count": len(affected_dgroup_ids),
                "affected_kegg_drug_count": len(affected_drug_ids),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(FAMILY_OUTPUT_DIR / "family_pharmacogenomics_summary.csv", index=False, encoding="utf-8-sig")

    print()
    print(summary_df.to_string(index=False))
    print()
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_pharmacogenomics.csv'}")
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_pharmacogenomics_summary.csv'}")


if __name__ == "__main__":
    main()
