"""kegg_relations.csv(엔트리 간 관계 edge list)를 그래프 관점에서 분석한다.

networkx 없이 pandas만으로:
  1. 타입-타입 간 관계 집계 (예: disease -> gene_entrez 몇 건)
  2. 노드별 연결 수(degree) 상위 허브 추출, 이름 라벨링
  3. 가장 연결이 많은 disease 노드 하나를 골라 ego network(이웃 서브그래프)를
     JSON으로 추출 -> 시각화 아티팩트에서 사용

출력 (results/tables/):
  kegg_graph_type_summary.csv
  kegg_graph_top_hubs.csv
  kegg_graph_ego_network.json
"""

import json
import sys
from pathlib import Path

import pandas as pd

from kegg_analysis_utils import build_gene_symbol_map

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RAW_KEGG_DIR = ROOT / "data" / "raw" / "kegg"
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "results" / "tables"

TOP_N_HUBS = 15

# 각 엔트리 타입별로 이름을 붙일 때 쓸 (파일, id컬럼, name컬럼)
LABEL_SOURCES = {
    "disease": ("kegg_disease.csv", "entry_id", "name"),
    "drug": ("kegg_drug.csv", "entry_id", "name"),
    "dgroup": ("kegg_dgroup.csv", "entry_id", "name"),
    "network": ("kegg_network.csv", "entry_id", "name"),
    "variant": ("kegg_variant.csv", "entry_id", "name"),
}


def load_labels():
    """entry_id -> name 매핑을 타입별로 만든다. ortholog/compound는 라벨 테이블이 없어 ID를 그대로 쓴다."""
    labels = {}
    for entity_type, (filename, id_col, name_col) in LABEL_SOURCES.items():
        df = pd.read_csv(PROCESSED_DIR / filename, dtype=str)
        labels[entity_type] = dict(zip(df[id_col], df[name_col].fillna("")))
    labels["gene_entrez"] = build_gene_symbol_map(RAW_KEGG_DIR)
    return labels


def label_of(labels, node_type, node_id):
    table = labels.get(node_type)
    if table and node_id in table and table[node_id]:
        name = table[node_id]
        return name if len(name) <= 60 else name[:57] + "..."
    return node_id


def build_type_summary(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["source_type", "target_type"])
        .size()
        .reset_index(name="edge_count")
        .sort_values("edge_count", ascending=False)
    )
    return summary


def build_degree_table(df: pd.DataFrame) -> pd.DataFrame:
    """source/target 양쪽에 등장한 횟수를 합쳐 (node_type, node_id) 별 degree를 구한다."""
    src = df[["source_type", "source_id"]].rename(columns={"source_type": "node_type", "source_id": "node_id"})
    dst = df[["target_type", "target_id"]].rename(columns={"target_type": "node_type", "target_id": "node_id"})
    all_nodes = pd.concat([src, dst], ignore_index=True)
    degree = all_nodes.groupby(["node_type", "node_id"]).size().reset_index(name="degree")
    return degree.sort_values("degree", ascending=False)


def build_top_hubs(degree: pd.DataFrame, labels: dict) -> pd.DataFrame:
    rows = []
    for node_type, group in degree.groupby("node_type"):
        top = group.sort_values("degree", ascending=False).head(TOP_N_HUBS)
        for _, row in top.iterrows():
            rows.append(
                {
                    "node_type": node_type,
                    "node_id": row["node_id"],
                    "label": label_of(labels, node_type, row["node_id"]),
                    "degree": row["degree"],
                }
            )
    return pd.DataFrame(rows)


def build_ego_network(df: pd.DataFrame, labels: dict, center_type: str, center_id: str, max_neighbors: int = 200):
    """center 노드의 1-hop 이웃 서브그래프를 뽑아 시각화용 JSON 구조로 만든다."""
    mask_out = (df["source_type"] == center_type) & (df["source_id"] == center_id)
    mask_in = (df["target_type"] == center_type) & (df["target_id"] == center_id)
    neighbor_edges = df[mask_out | mask_in].copy()

    if len(neighbor_edges) > max_neighbors:
        neighbor_edges = neighbor_edges.head(max_neighbors)

    nodes = {}
    center_key = f"{center_type}:{center_id}"
    nodes[center_key] = {
        "id": center_key,
        "type": center_type,
        "label": label_of(labels, center_type, center_id),
        "is_center": True,
    }

    edges = []
    for _, row in neighbor_edges.iterrows():
        if row["source_type"] == center_type and row["source_id"] == center_id:
            other_type, other_id = row["target_type"], row["target_id"]
        else:
            other_type, other_id = row["source_type"], row["source_id"]

        other_key = f"{other_type}:{other_id}"
        if other_key not in nodes:
            nodes[other_key] = {
                "id": other_key,
                "type": other_type,
                "label": label_of(labels, other_type, other_id),
                "is_center": False,
            }
        edges.append({"source": center_key, "target": other_key, "field": row["field"]})

    return {"center": center_key, "nodes": list(nodes.values()), "edges": edges}


def main():
    df = pd.read_csv(PROCESSED_DIR / "kegg_relations.csv", dtype=str)
    labels = load_labels()

    type_summary = build_type_summary(df)
    type_summary.to_csv(OUTPUT_DIR / "kegg_graph_type_summary.csv", index=False, encoding="utf-8-sig")

    degree = build_degree_table(df)
    top_hubs = build_top_hubs(degree, labels)
    top_hubs.to_csv(OUTPUT_DIR / "kegg_graph_top_hubs.csv", index=False, encoding="utf-8-sig")

    disease_degree = degree[degree["node_type"] == "disease"].sort_values("degree", ascending=False)
    center_id = disease_degree.iloc[0]["node_id"]
    ego = build_ego_network(df, labels, "disease", center_id)

    with open(OUTPUT_DIR / "kegg_graph_ego_network.json", "w", encoding="utf-8") as f:
        json.dump(ego, f, ensure_ascii=False, indent=2)

    print("타입 간 관계 요약 (상위 5):")
    print(type_summary.head(5).to_string(index=False))
    print()
    print("가장 연결이 많은 disease 허브:", center_id, "-", label_of(labels, "disease", center_id))
    print(f"  degree={int(disease_degree.iloc[0]['degree'])}, ego network 노드 {len(ego['nodes'])}개")
    print()
    print(f"-> {OUTPUT_DIR / 'kegg_graph_type_summary.csv'}")
    print(f"-> {OUTPUT_DIR / 'kegg_graph_top_hubs.csv'}")
    print(f"-> {OUTPUT_DIR / 'kegg_graph_ego_network.json'}")


if __name__ == "__main__":
    main()
