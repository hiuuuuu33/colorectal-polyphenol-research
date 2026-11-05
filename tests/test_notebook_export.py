"""Check that HTML export preserves R syntax and excludes formula processing."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from bs4 import BeautifulSoup
import nbformat


EXPORTER = Path(__file__).resolve().parents[1] / "scripts/render_notebook.py"


class NotebookExportTest(unittest.TestCase):
    def test_r_code_is_literal_and_notebook_is_unchanged(self):
        sources = [
            'na <- de$Adenoma_vs_Normal\n'
            'up <- intersect(na$gene[na$direction == "Up"], ac$gene)\n'
            'logFC_NA <- na$logFC[match(genes, na$gene)]',
            'label <- "$value_with_underscores$"\n'
            'path <- "C:\\\\data\\\\file.csv"\n'
            'if (x < 1 && y > 0) {\n  print("<tag> & text")\n}',
        ]
        notebook = nbformat.v4.new_notebook(
            metadata={"language_info": {"name": "R", "pygments_lexer": "r"}},
            cells=[nbformat.v4.new_code_cell(source, execution_count=i + 1,
                                             id=f"syntax-{i}")
                   for i, source in enumerate(sources)],
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "report").mkdir()
            (root / "report/notebook.css").write_text("")
            path = root / "Analysis.ipynb"
            nbformat.write(notebook, path)
            original = path.read_bytes()
            subprocess.run([sys.executable, str(EXPORTER), "--root", str(root),
                            "--export-only"], check=True, capture_output=True, text=True)
            soup = BeautifulSoup((root / "docs/index.html").read_text(), "html.parser")
            self.assertFalse(soup.find_all("script"))
            cells = soup.select(".code_cell")
            self.assertEqual(len(cells), len(sources))
            for i, (cell, source) in enumerate(zip(cells, sources)):
                self.assertEqual(cell["id"], f"cell-id=syntax-{i}")
                area = cell.select_one(".input_area.tex2jax_ignore.mathjax_ignore")
                self.assertIsNotNone(area)
                code = area.select_one("pre > code")
                self.assertIsNotNone(code)
                self.assertEqual(code.get_text(), source + "\n")
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
