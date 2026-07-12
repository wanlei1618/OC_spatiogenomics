from pathlib import Path
from zipfile import ZipFile
from lxml import etree

docx_path = Path(r"D:\Downloads\CNV_expression_joint_analysis_integrated_oc (1).docx")
out_path = Path(r"C:\Users\chenfy12\Documents\Codex\2026-07-05\list-files-home-shpc-006-oc\work\CNV_expression_joint_analysis_integrated_oc.txt")
ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}

def paragraph_text(p):
    parts = []
    for node in p.iter():
        if node.tag == f"{{{ns['w']}}}t":
            parts.append(node.text or "")
        elif node.tag == f"{{{ns['w']}}}tab":
            parts.append("\t")
        elif node.tag == f"{{{ns['w']}}}br":
            parts.append("\n")
    return "".join(parts).strip()

with ZipFile(docx_path) as zf:
    root = etree.fromstring(zf.read("word/document.xml"))

lines = []
for child in root.xpath("//w:body/*", namespaces=ns):
    if child.tag == f"{{{ns['w']}}}p":
        text = paragraph_text(child)
        if text:
            lines.append(text)
    elif child.tag == f"{{{ns['w']}}}tbl":
        rows = []
        for tr in child.xpath(".//w:tr", namespaces=ns):
            cells = []
            for tc in tr.xpath("./w:tc", namespaces=ns):
                ps = [paragraph_text(p) for p in tc.xpath(".//w:p", namespaces=ns)]
                cells.append(" ".join([p for p in ps if p]))
            rows.append("\t".join(cells))
        if rows:
            lines.append("\n".join(rows))

out_path.write_text("\n\n".join(lines), encoding="utf-8")
print(out_path)
