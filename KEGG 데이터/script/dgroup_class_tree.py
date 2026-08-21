"""dgroup.txt의 CLASS 필드를 파싱해 약물 분류 트리(DAG)를 만든다.

CLASS 필드는 이런 식으로 되어 있다 (예: DG00004 Miconazole):
    CLASS       Antifungal
                 DG01883  Imidazole antifungal
                Metabolizing enzyme inhibitor
                 DG01643  CYP2C9 inhibitor
                 DG02852  CYP3A/CYP3A4 inhibitor

들여쓰기가 없는 줄(예: "Antifungal")은 최상위 텍스트 카테고리(고정 어휘,
KEGG가 애초에 DG 번호를 안 붙인 루트 노드)이고, "DG#####  이름" 형식인 줄은
그 카테고리 아래의 특정 dgroup을 가리키는 부모-자식 엣지다. 한 dgroup이 여러
카테고리에 동시에 속할 수 있어(예: 항진균제이면서 대사효소억제제이기도 함)
엄밀한 트리가 아니라 다중부모 DAG다.

출력 (KEGG 데이터/output/):
  dgroup_class_edges.csv  - (child_dgroup_id, parent_type, parent_id, parent_label) 엣지 목록
  dgroup_class_rollup.csv - 루트 카테고리별 하위 dgroup 수·고유 약물 수(재귀 집계)
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kegg_flatfile_parser import parse_entries, entry_id_and_type, field_lines  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_DIR = SCRIPT_DIR.parent / "output"

DG_LINE_RE = re.compile(r"^(DG\d{5})\s+(.*)$")


def parse_class_edges(entry_id: str, class_lines: list):
    """CLASS 필드 줄들을 (parent_type, parent_id, parent_label) 엣지로 변환한다.
    루트 텍스트 카테고리 뒤에 DG 줄이 하나도 안 따라오면 그 루트 자체가 직접 부모가 된다."""
    edges = []
    current_root = None
    root_has_dg_child = {}

    for line in class_lines:
        m = DG_LINE_RE.match(line)
        if m:
            dg_id, dg_name = m.group(1), m.group(2)
            edges.append(("dgroup", dg_id, dg_name))
            if current_root is not None:
                root_has_dg_child[current_root] = True
        else:
            current_root = line
            root_has_dg_child.setdefault(current_root, False)

    for root, has_child in root_has_dg_child.items():
        if not has_child:
            edges.append(("root", root, root))

    return [(entry_id, ptype, pid, plabel) for ptype, pid, plabel in edges]


def main():
    entries = parse_entries(DATA_DIR / "dgroup.txt")

    dgroup_name = {}
    dgroup_member_ids = {}
    all_edges = []

    for entry in entries:
        entry_id, _ = entry_id_and_type(entry)
        if not entry_id:
            continue
        name_lines = entry.get("NAME", [])
        dgroup_name[entry_id] = name_lines[0].strip() if name_lines else entry_id

        member_ids = re.findall(r"\bD\d{5}\b", " ".join(entry.get("MEMBER", [])))
        dgroup_member_ids[entry_id] = set(member_ids)

        all_edges.extend(parse_class_edges(entry_id, field_lines(entry, "CLASS")))

    edges_df = pd.DataFrame(all_edges, columns=["child_dgroup_id", "parent_type", "parent_id", "parent_label"])
    edges_df["child_name"] = edges_df["child_dgroup_id"].map(dgroup_name)
    edges_df.to_csv(OUTPUT_DIR / "dgroup_class_edges.csv", index=False, encoding="utf-8-sig")
    print(f"CLASS 엣지 {len(edges_df)}건 ({edges_df['child_dgroup_id'].nunique()}개 dgroup에 CLASS 정보 있음)")

    # 자식 -> 부모 인접 리스트 (부모가 root거나 dgroup이거나 상관없이)
    children_of = defaultdict(set)  # parent_key -> {child_dgroup_id}
    for child, ptype, pid, plabel in all_edges:
        parent_key = ("root", plabel) if ptype == "root" else ("dgroup", pid)
        children_of[parent_key].add(child)

    def collect_descendants(parent_key, seen=None):
        """parent 아래 모든 dgroup(자기 자신 포함)을 재귀로 모은다. 순환 방지용 seen 집합 사용."""
        if seen is None:
            seen = set()
        for child in children_of.get(parent_key, ()):
            if child in seen:
                continue
            seen.add(child)
            collect_descendants(("dgroup", child), seen)
        return seen

    root_labels = sorted({plabel for _, ptype, _, plabel in all_edges if ptype == "root"})
    rollup_rows = []
    for root in root_labels:
        descendants = collect_descendants(("root", root))
        member_union = set()
        for d in descendants:
            member_union |= dgroup_member_ids.get(d, set())
        rollup_rows.append(
            {
                "root_category": root,
                "descendant_dgroup_count": len(descendants),
                "unique_drug_count": len(member_union),
            }
        )

    rollup_df = pd.DataFrame(rollup_rows).sort_values("unique_drug_count", ascending=False)
    rollup_df.to_csv(OUTPUT_DIR / "dgroup_class_rollup.csv", index=False, encoding="utf-8-sig")

    print()
    print(rollup_df.head(20).to_string(index=False))
    print()
    print(f"-> {OUTPUT_DIR / 'dgroup_class_edges.csv'}")
    print(f"-> {OUTPUT_DIR / 'dgroup_class_rollup.csv'}")


if __name__ == "__main__":
    main()
