"""Verify the extracted release with the Python standard library."""
from pathlib import Path
import csv
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET

root = Path(__file__).resolve().parent
manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
for row in manifest:
    p = root / row["file"]
    assert p.is_file(), row["file"]
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row["sha256"], row["file"]
for p in root.rglob("*.xlsx"):
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        if "docProps/core.xml" in z.namelist():
            core = ET.fromstring(z.read("docProps/core.xml"))
            assert not any(x.tag.rsplit("}", 1)[-1] in
                           {"creator", "lastModifiedBy", "created", "modified"} for x in core)
with (root / "original_analysis/data/processed/singleton_analysis_with_tail_phenotypes.csv").open(
        encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 9384
assert len({r["dam_id"] for r in rows}) == 4624
assert sum(int(r["lower_tail"]) for r in rows) == 948
assert sum(int(r["upper_tail"]) for r in rows) == 955
assert len(list((root / "figures").glob("*.png"))) == 4
with zipfile.ZipFile(root / "tables/Supplementary_Tables.xlsx") as z:
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    book = ET.fromstring(z.read("xl/workbook.xml"))
    assert [x.attrib["name"] for x in book.findall("s:sheets/s:sheet", ns)] == [
        "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
print(f"Verified {len(manifest)} file hashes, metadata, cohort totals, four figures and seven tables.")
