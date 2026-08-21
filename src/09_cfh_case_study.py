"""가족이 실제로 보유한 CFH 병원성 변이(Phase 3)를 Phase 4의 질병 중심 분석 엔진에 그대로
대입해, "이 가족에게 구체적으로 어떤 의미인지"를 하나의 프로파일로 만든다.

CFH(Complement Factor H) 변이는 disease.txt 기준 5개 질병과 연결돼 있다:
  - Alternative complement pathway component defects (H00104)
  - Age-related macular degeneration (H00821)
  - Atypical hemolytic uremic syndrome (H01434)
  - Basal laminar drusen (H02108)
  - C3 glomerulopathy (H02579)

이 5개 질병에 한해 disease_centered_analysis.py와 같은 로직으로 약물 재창출 후보와
동반질환 쌍을 다시 계산한다 — 전역 top-300 컷오프 없이 이 5개 질병에 한정해서 전부 본다.

출력 (results/tables/):
  cfh_family_disease_profile.csv - CFH 연결 질병별 요약 + 재창출 후보 + 동반질환
"""

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from kegg_analysis_utils import build_gene_symbol_map, build_index, load_name_maps, short

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RAW_KEGG_DIR = ROOT / "data" / "raw" / "kegg"
PROCESSED_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "results" / "tables"

CFH_DISEASE_IDS = ["H00104", "H00821", "H01434", "H02108", "H02579"]


def main():
    rel = pd.read_csv(PROCESSED_DIR / "kegg_relations.csv", dtype=str)
    names = load_name_maps(PROCESSED_DIR)
    gene_symbol = build_gene_symbol_map(RAW_KEGG_DIR)

    disease_genes_direct = build_index(rel, "disease", "gene_entrez")
    drug_gene = build_index(rel, "drug", "gene_entrez")
    gene_to_drugs = defaultdict(set)
    for drug_id, gids in drug_gene.items():
        for gid in gids:
            gene_to_drugs[gid].add(drug_id)
    disease_drugs_known = build_index(rel, "disease", "drug")
    network_diseases = build_index(rel, "network", "disease")
    disease_networks = defaultdict(set)
    for nid, dids in network_diseases.items():
        for did in dids:
            disease_networks[did].add(nid)

    gene_diseases = defaultdict(set)
    for did, gids in disease_genes_direct.items():
        for gid in gids:
            gene_diseases[gid].add(did)

    rows = []
    for did in CFH_DISEASE_IDS:
        direct_genes = disease_genes_direct.get(did, set())
        gene_syms = sorted({gene_symbol.get(g, g) for g in direct_genes})

        candidate_genes = defaultdict(set)
        for gid in direct_genes:
            for drug_id in gene_to_drugs.get(gid, ()):
                candidate_genes[drug_id].add(gid)
        known = disease_drugs_known.get(did, set())
        novel = {d: g for d, g in candidate_genes.items() if d not in known}
        top_candidates = sorted(novel.items(), key=lambda kv: -len(kv[1]))[:8]
        candidate_text = "; ".join(
            f"{names['drug_name'].get(d, d).split(';')[0].strip()} (via {'/'.join(gene_symbol.get(g, g) for g in g_)})"
            for d, g_ in top_candidates
        )

        # 이 5개 질병 간 + 이 5개 질병과 다른 모든 질병 간 동반질환(공유 유전자) 쌍
        comorbid = []
        for other_did in gene_diseases_of(direct_genes, gene_diseases):
            if other_did == did:
                continue
            other_genes = disease_genes_direct.get(other_did, set())
            shared = direct_genes & other_genes
            if shared:
                comorbid.append((other_did, len(shared)))
        comorbid.sort(key=lambda x: -x[1])
        comorbid_text = "; ".join(
            f"{short(names['disease_name'].get(o, o), 40)} (공유 유전자 {c}개)" for o, c in comorbid[:6]
        )

        rows.append(
            {
                "disease_id": did,
                "disease_name": names["disease_name"].get(did, did),
                "cfh_pathway_genes": ", ".join(gene_syms),
                "drug_repurposing_candidates": candidate_text,
                "comorbid_diseases_via_shared_gene": comorbid_text,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "cfh_family_disease_profile.csv", index=False, encoding="utf-8-sig")

    pd.set_option("display.max_colwidth", 60)
    print(df.to_string(index=False))
    print()
    print(f"-> {RESULTS_DIR / 'cfh_family_disease_profile.csv'}")


def gene_diseases_of(gene_ids, gene_diseases_index):
    result = set()
    for g in gene_ids:
        result |= gene_diseases_index.get(g, set())
    return result


if __name__ == "__main__":
    main()
