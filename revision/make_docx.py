#!/usr/bin/env python3
"""make_docx.py — build PLOS-style docx files from the Markdown deliverables.

Usage: python3 -I make_docx.py <kind> <in.md> <out.docx> [--linenum]
  kind = manuscript | letter | appendix
Requires pandoc and python-docx.
"""
import sys, os, re, subprocess, zipfile, shutil, tempfile
from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

kind, src, out = sys.argv[1:4]
LINENUM = "--linenum" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))

md = open(src, encoding="utf-8").read()
md = re.sub(r"<!--.*?-->", "", md, flags=re.S)
md = md.split("\n# 先生への確認事項")[0]
if kind == "appendix":
    md = re.sub(r"\s*⟨図 (\S+?) を挿入⟩", r"\n\n[Figure file \1 to be inserted here]\n", md)
md = re.sub(r"\s*⟨[^⟩]*⟩", "", md)
if kind == "appendix":
    md = re.sub(r"^# S1 Appendix v4.*$", "# S1 Appendix. Supplementary methods and results", md, flags=re.M)
if kind == "letter":
    md = re.sub(r"^# Response to Reviewers.*$", "# Response to Reviewers — PDIG-D-26-01140", md, flags=re.M)
    md = re.sub(r"^\*\*Status \(.*?\n", "", md, flags=re.M)
    md = re.sub(r"^\((i|ii|iii|iv|v)\) ", lambda m: "(" + m.group(1) + ") ", md, flags=re.M)
tmp = tempfile.mkdtemp()
clean = os.path.join(tmp, "clean.md"); open(clean, "w", encoding="utf-8").write(md)

ref_dir = os.path.join(tmp, "ref"); os.makedirs(ref_dir)
subprocess.run(["pandoc", "-o", os.path.join(tmp, "ref.docx"), "--print-default-data-file", "reference.docx"], check=True)
with zipfile.ZipFile(os.path.join(tmp, "ref.docx")) as z: z.extractall(ref_dir)
sty = open(os.path.join(ref_dir, "word/styles.xml"), encoding="utf-8").read()
base = {"manuscript": 24, "letter": 22, "appendix": 20}[kind]
line = {"manuscript": 480, "letter": 276, "appendix": 276}[kind]
font = '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" w:cs="Times New Roman" />'
sty = re.sub(r'<w:rFonts w:asciiTheme="minorHAnsi"[^>]*/>', font, sty, count=1)
sty = sty.replace('<w:sz w:val="24" />\n        <w:szCs w:val="24" />', f'<w:sz w:val="{base}" />\n        <w:szCs w:val="{base}" />', 1)
sty = sty.replace('<w:pPr>\n        <w:spacing w:after="200" />\n      </w:pPr>', f'<w:pPr>\n        <w:spacing w:after="120" w:line="{line}" w:lineRule="auto" />\n      </w:pPr>')
sty = re.sub(r'<w:rFonts w:asciiTheme="majorHAnsi"[^>]*/>', font, sty)
sty = re.sub(r'<w:color w:val="[0-9A-F]{6}" w:themeColor="accent1" />', '<w:color w:val="000000" />', sty)
sty = re.sub(r'(<w:style w:type="paragraph" w:styleId="Heading[1-9]">.*?)<w:sz w:val="\d+" />', lambda m: m.group(1) + f'<w:sz w:val="{base}" />', sty, flags=re.S)
sty = re.sub(r'(<w:style w:type="paragraph" w:styleId="Heading[1-9]">.*?<w:rPr>.*?)<w:i />\s*<w:iCs />', r'\1', sty, flags=re.S)
sty = re.sub(r'(<w:style w:type="paragraph" w:styleId="Heading[1-9]">.*?)<w:color w:val="[0-9A-F]{6}" w:themeColor="\w+"( w:themeShade="\w+")? />', r'\1<w:color w:val="000000" />', sty, flags=re.S)
if kind != "manuscript":
    sty = sty.replace('<w:name w:val="Compact" />\n    <w:basedOn w:val="BodyText" />\n    <w:qFormat />\n    <w:pPr>\n      <w:spacing w:before="36" w:after="36" />\n    </w:pPr>',
                      f'<w:name w:val="Compact" />\n    <w:basedOn w:val="BodyText" />\n    <w:qFormat />\n    <w:pPr>\n      <w:spacing w:before="20" w:after="20" w:line="240" w:lineRule="auto" />\n    </w:pPr>\n    <w:rPr><w:sz w:val="{16 if kind=="appendix" else 18}" /><w:szCs w:val="{16 if kind=="appendix" else 18}" /></w:rPr>')
else:
    sty = sty.replace('<w:spacing w:before="36" w:after="36" />', '<w:spacing w:before="20" w:after="20" w:line="240" w:lineRule="auto" />')
open(os.path.join(ref_dir, "word/styles.xml"), "w", encoding="utf-8").write(sty)
if kind == "appendix":
    d = open(os.path.join(ref_dir, "word/document.xml"), encoding="utf-8").read()
    sect = '<w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape" /><w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" w:header="708" w:footer="708" w:gutter="0" /></w:sectPr>'
    d = re.sub(r'<w:sectPr\b.*?</w:sectPr>|<w:sectPr\b[^>]*/>', '', d, flags=re.S)
    d = d.replace('</w:body>', sect + '</w:body>')
    open(os.path.join(ref_dir, "word/document.xml"), "w", encoding="utf-8").write(d)
ref = os.path.join(tmp, "reference.docx")
with zipfile.ZipFile(ref, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(ref_dir):
        for f in files:
            p = os.path.join(root, f); z.write(p, os.path.relpath(p, ref_dir))

subprocess.run(["pandoc", clean, "--reference-doc", ref, "--resource-path", os.path.dirname(os.path.abspath(src)), "-o", out], check=True)

doc = Document(out)
if kind == "manuscript":  # numbered/bulleted list items (Compact style) double-spaced like the body text; tables keep single spacing
    for p in doc.paragraphs:
        if p._p.pPr is not None and p._p.pPr.find(qn("w:numPr")) is not None:
            sp = p._p.pPr.find(qn("w:spacing"))
            if sp is None:
                sp = OxmlElement("w:spacing"); p._p.pPr.append(sp)
            sp.set(qn("w:line"), "480"); sp.set(qn("w:lineRule"), "auto"); sp.set(qn("w:before"), "0"); sp.set(qn("w:after"), "120")
if LINENUM:
    for sect in doc.sections:
        sp = sect._sectPr
        ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "360")
        sp.insert(0, ln)
page_w = int((Cm(29.7 - 4) if kind == "appendix" else Cm(21.0 - 5.08)).emu * 0.96)
def set_borders(tbl):
    tblPr = tbl._tbl.tblPr
    b = OxmlElement("w:tblBorders")
    for edge in ("top", "bottom", "insideH"):
        e = OxmlElement(f"w:{edge}"); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "4"); e.set(qn("w:color"), "000000"); b.append(e)
    tblPr.append(b)
for tbl in doc.tables:
    ncol = len(tbl.columns)
    lens = [1] * ncol
    for row in tbl.rows:
        for i, c in enumerate(row.cells[:ncol]):
            lens[i] = max(lens[i], min(max(len(c.text), 14), 60))
    tot = sum(lens)
    widths = [max(int(page_w * l / tot), Cm(1.2).emu) for l in lens]
    scale = page_w / sum(widths); widths = [int(w * scale) for w in widths]
    tbl.autofit = False
    tblPr = tbl._tbl.tblPr
    for old in tblPr.findall(qn("w:tblW")) + tblPr.findall(qn("w:tblLayout")): tblPr.remove(old)
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblPr.append(lay)
    tw = OxmlElement("w:tblW"); tw.set(qn("w:w"), str(int(page_w / 635))); tw.set(qn("w:type"), "dxa"); tblPr.append(tw)
    grid = tbl._tbl.find(qn("w:tblGrid"))
    if grid is not None:
        for i, gc in enumerate(grid.findall(qn("w:gridCol"))[:ncol]):
            gc.set(qn("w:w"), str(int(widths[i] / 635)))
    for row in tbl.rows:
        for i, c in enumerate(row.cells[:ncol]):
            c.width = widths[i]
            for p in c.paragraphs:
                for r in p.runs:
                    r.font.name = "Times New Roman"
                    r.font.size = Pt(8 if (kind == "appendix" or ncol > 6) else 10)
    set_borders(tbl)
for p in doc.paragraphs:
    for h in p._p.findall(".//" + qn("w:hyperlink")):
        for r in h.findall(".//" + qn("w:rPr")):
            for c in r.findall(qn("w:color")): c.set(qn("w:val"), "000000")
doc.save(out)
shutil.rmtree(tmp)
print("wrote", out)
