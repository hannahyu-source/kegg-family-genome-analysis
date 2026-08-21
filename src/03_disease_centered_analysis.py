"""질병 중심 분석 — kegg_relations.csv를 이용한 두 가지 분석.

1. 약물 재창출(drug repurposing) 후보 탐색
   질병 -> (직접 연결된 유전자) + (연결된 network의 유전자) = "병리 경로 유전자 집합"
   그 유전자를 표적(TARGET)으로 하는 약물 중, 이 질병에 아직 공식적으로 연결(disease->drug)
   되지 않은 약물을 "후보"로 뽑는다. 공유 유전자가 많을수록 더 강한 후보로 본다.

   주의: network.txt는 개별 network 엔트리(N#####)의 DISEASE 필드에 관련 질병을 직접
   명시한다. disease.txt의 NETWORK 필드는 반대로 nt0xxxx CLASS 카테고리 코드를 적는데,
   이건 N#####과 다른 네임스페이스라 엔트리 단위로 직접 연결되지 않는다. 그래서 "질병에
   연결된 network"는 network.txt DISEASE 필드(kegg_relations의 network->disease 엣지)
   방향만 신뢰할 수 있는 소스로 쓴다.

2. 동반질환(comorbidity) 패턴 분석
   같은 network 또는 같은 유전자를 공유하는 질병 쌍을 찾아 "기전적으로 연결된" 후보로 본다.
   (진짜 임상적 동반이환 여부가 아니라, 분자 기전을 공유한다는 뜻 — 보고서에도 이 점을 명시)

출력 (results/tables/):
  disease_pathway_genes.csv       - 질병별 관련 network 수 / 병리경로 유전자 수
  drug_repurposing_candidates.csv - 질병별 상위 15개 재창출 후보 약물
  disease_comorbidity_pairs.csv   - 공유 network/유전자 기준 질병 쌍 랭킹
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
OUTPUT_DIR = ROOT / "results" / "tables"

TOP_CANDIDATES_PER_DISEASE = 15
TOP_COMORBIDITY_PAIRS = 300


def build_repurposing_candidates(rel, names, gene_symbol):
    network_diseases = build_index(rel, "network", "disease")
    disease_networks = defaultdict(set)
    for nid, dids in network_diseases.items():
        for did in dids:
            disease_networks[did].add(nid)

    disease_genes_direct = build_index(rel, "disease", "gene_entrez")
    network_genes = build_index(rel, "network", "gene_entrez")
    drug_gene = build_index(rel, "drug", "gene_entrez")
    gene_to_drugs = defaultdict(set)
    for drug_id, gids in drug_gene.items():
        for gid in gids:
            gene_to_drugs[gid].add(drug_id)

    disease_drugs_known = build_index(rel, "disease", "drug")

    summary_rows = []
    candidate_rows = []

    for did in names["disease_ids"]:
        networks = disease_networks.get(did, set())
        direct_genes = disease_genes_direct.get(did, set())
        pathway_genes = set(direct_genes)
        for nid in networks:
            pathway_genes |= network_genes.get(nid, set())

        summary_rows.append(
            {
                "disease_id": did,
                "disease_name": short(names["disease_name"].get(did, did)),
                "linked_network_count": len(networks),
                "direct_gene_count": len(direct_genes),
                "pathway_gene_count": len(pathway_genes),
            }
        )

        if not pathway_genes:
            continue

        candidate_genes = defaultdict(set)  # drug_id -> matched gene ids
        for gid in pathway_genes:
            for drug_id in gene_to_drugs.get(gid, ()):
                candidate_genes[drug_id].add(gid)

        known = disease_drugs_known.get(did, set())
        novel = {d: g for d, g in candidate_genes.items() if d not in known}
        top = sorted(novel.items(), key=lambda kv: -len(kv[1]))[:TOP_CANDIDATES_PER_DISEASE]

        for drug_id, gids in top:
            gene_syms = sorted({gene_symbol.get(g, g) for g in gids})
            candidate_rows.append(
                {
                    "disease_id": did,
                    "disease_name": short(names["disease_name"].get(did, did)),
                    "candidate_drug_id": drug_id,
                    "candidate_drug_name": short(names["drug_name"].get(drug_id, drug_id)),
                    "candidate_dgroup_id": names["drug_dgroup"].get(drug_id, ""),
                    "shared_gene_count": len(gids),
                    "shared_gene_symbols": "; ".join(gene_syms),
                }
            )

    summary_df = pd.DataFrame(summary_rows).sort_values("pathway_gene_count", ascending=False)
    candidates_df = pd.DataFrame(candidate_rows).sort_values(
        ["disease_id", "shared_gene_count"], ascending=[True, False]
    )
    return summary_df, candidates_df


def build_comorbidity_pairs(rel, names):
    network_diseases = build_index(rel, "network", "disease")
    disease_networks = defaultdict(set)
    for nid, dids in network_diseases.items():
        for did in dids:
            disease_networks[did].add(nid)

    disease_genes_direct = build_index(rel, "disease", "gene_entrez")
    gene_diseases = defaultdict(set)
    for did, gids in disease_genes_direct.items():
        for gid in gids:
            gene_diseases[gid].add(did)

    pair_networks = defaultdict(set)
    for nid, dids in network_diseases.items():
        if len(dids) < 2:
            continue
        dl = sorted(dids)
        for i in range(len(dl)):
            for j in range(i + 1, len(dl)):
                pair_networks[(dl[i], dl[j])].add(nid)

    pair_genes = defaultdict(set)
    for gid, dids in gene_diseases.items():
        if len(dids) < 2:
            continue
        dl = sorted(dids)
        for i in range(len(dl)):
            for j in range(i + 1, len(dl)):
                pair_genes[(dl[i], dl[j])].add(gid)

    all_pairs = set(pair_networks) | set(pair_genes)
    rows = []
    for d1, d2 in all_pairs:
        shared_n = pair_networks.get((d1, d2), set())
        shared_g = pair_genes.get((d1, d2), set())
        rows.append(
            {
                "disease_a_id": d1,
                "disease_a_name": short(names["disease_name"].get(d1, d1)),
                "disease_b_id": d2,
                "disease_b_name": short(names["disease_name"].get(d2, d2)),
                "shared_network_count": len(shared_n),
                "shared_gene_count": len(shared_g),
                "shared_networks": "; ".join(short(names["network_name"].get(n, n), 40) for n in sorted(shared_n)),
                "score": len(shared_n) + len(shared_g),
            }
        )

    df = pd.DataFrame(rows).sort_values("score", ascending=False).drop(columns="score")
    return df.head(TOP_COMORBIDITY_PAIRS)


def main():
    rel = pd.read_csv(PROCESSED_DIR / "kegg_relations.csv", dtype=str)
    names = load_name_maps(PROCESSED_DIR)
    gene_symbol = build_gene_symbol_map(RAW_KEGG_DIR)

    print("1) 약물 재창출 후보 탐색 중...")
    summary_df, candidates_df = build_repurposing_candidates(rel, names, gene_symbol)
    summary_df.to_csv(OUTPUT_DIR / "disease_pathway_genes.csv", index=False, encoding="utf-8-sig")
    candidates_df.to_csv(OUTPUT_DIR / "drug_repurposing_candidates.csv", index=False, encoding="utf-8-sig")
    print(f"   질병 {len(summary_df)}개, 재창출 후보 {len(candidates_df)}건")

    print("2) 동반질환 패턴(공유 network/유전자) 분석 중...")
    comorbidity_df = build_comorbidity_pairs(rel, names)
    comorbidity_df.to_csv(OUTPUT_DIR / "disease_comorbidity_pairs.csv", index=False, encoding="utf-8-sig")
    print(f"   상위 {len(comorbidity_df)}개 질병 쌍 저장")

    print()
    print("=== 예시: Hepatocellular carcinoma(H00048) 재창출 후보 상위 5 ===")
    ex = candidates_df[candidates_df.disease_id == "H00048"].head(5)
    print(ex[["candidate_drug_name", "shared_gene_count", "shared_gene_symbols"]].to_string(index=False))

    print()
    print("=== 동반질환 패턴 상위 5쌍 ===")
    print(comorbidity_df.head(5)[["disease_a_name", "disease_b_name", "shared_network_count", "shared_gene_count"]].to_string(index=False))

    print()
    print(f"-> {OUTPUT_DIR / 'disease_pathway_genes.csv'}")
    print(f"-> {OUTPUT_DIR / 'drug_repurposing_candidates.csv'}")
    print(f"-> {OUTPUT_DIR / 'disease_comorbidity_pairs.csv'}")


if __name__ == "__main__":
    main()
