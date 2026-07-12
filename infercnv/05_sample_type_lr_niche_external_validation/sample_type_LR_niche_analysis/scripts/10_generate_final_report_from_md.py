import html
import re
import zipfile
from pathlib import Path

import markdown


ROOT = Path(r"D:\OC_spatiogenomics\infercnv\sample_type_LR_niche_analysis")
TABLE_DIR = ROOT / "tables"
MD = ROOT / "sample_type_LR_niche_analysis_report.md"
HTML = ROOT / "sample_type_LR_niche_analysis_report.html"
DOCX = ROOT / "sample_type_LR_niche_analysis_report.docx"


def markdown_to_html():
    text = MD.read_text(encoding="utf-8")
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    page = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>sample_type LR niche analysis</title>
  <style>
    body {{ font-family: Arial, sans-serif; max-width: 1100px; margin: 32px auto; line-height: 1.45; }}
    table {{ border-collapse: collapse; font-size: 12px; }}
    th, td {{ border: 1px solid #ccc; padding: 4px 6px; }}
    code {{ background: #f4f4f4; padding: 1px 3px; }}
    pre {{ background: #f7f7f7; padding: 12px; overflow-x: auto; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""
    HTML.write_text(page, encoding="utf-8")


def paragraph_xml(text, style=None):
    text = html.escape(text)
    p_style = ""
    if style:
        p_style = f"<w:pPr><w:pStyle w:val=\"{style}\"/></w:pPr>"
    return f"<w:p>{p_style}<w:r><w:t xml:space=\"preserve\">{text}</w:t></w:r></w:p>"


def markdown_to_docx():
    lines = MD.read_text(encoding="utf-8").splitlines()
    paragraphs = []
    for line in lines:
        raw = line.rstrip()
        if not raw:
            paragraphs.append("<w:p/>")
            continue
        if raw.startswith("# "):
            paragraphs.append(paragraph_xml(raw[2:], "Heading1"))
        elif raw.startswith("## "):
            paragraphs.append(paragraph_xml(raw[3:], "Heading2"))
        elif raw.startswith("### "):
            paragraphs.append(paragraph_xml(raw[4:], "Heading3"))
        else:
            text = re.sub(r"`([^`]*)`", r"\1", raw)
            paragraphs.append(paragraph_xml(text))

    document_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    {''.join(paragraphs)}
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"/>
    </w:sectPr>
  </w:body>
</w:document>
"""
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
"""
    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:styleId="Normal"><w:name w:val="Normal"/></w:style>
  <w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:outlineLvl w:val="0"/></w:pPr><w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/><w:pPr><w:outlineLvl w:val="1"/></w:pPr><w:rPr><w:b/><w:sz w:val="26"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading3"><w:name w:val="heading 3"/><w:basedOn w:val="Normal"/><w:pPr><w:outlineLvl w:val="2"/></w:pPr><w:rPr><w:b/><w:sz w:val="22"/></w:rPr></w:style>
</w:styles>
"""
    with zipfile.ZipFile(DOCX, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("word/document.xml", document_xml)
        z.writestr("word/styles.xml", styles)


def update_status():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    status = TABLE_DIR / "final_report_render_status.csv"
    status.write_text(
        "output,status,note\n"
        "html,completed_python_markdown_fallback,pandoc unavailable in R environment\n"
        "docx,completed_minimal_ooxml_fallback,pandoc and python-docx unavailable; generated valid minimal OOXML docx\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    markdown_to_html()
    markdown_to_docx()
    update_status()
