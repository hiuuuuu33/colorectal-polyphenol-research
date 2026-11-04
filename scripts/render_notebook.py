"""Execute the R notebook and export its cells and results as an offline HTML page."""

import argparse
import base64
import io
import os
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
from PIL import Image


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--kernel", default="crc-polyphenol-r")
    parser.add_argument("--export-only", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    notebook_path = root / "Analysis.ipynb"
    notebook = nbformat.read(notebook_path, as_version=4)
    nbformat.validate(notebook)
    if not args.export_only:
        os.environ["PYTHON"] = sys.executable
        os.environ["RESEARCH_REPO_ROOT"] = str(root)
        for cell in notebook.cells:
            if cell.cell_type == "code":
                cell.outputs = []
                cell.execution_count = None
        client = NotebookClient(notebook, timeout=600, kernel_name=args.kernel,
                                resources={"metadata": {"path": str(root)}})
        client.execute()
        nbformat.validate(notebook)
    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    if not code_cells or any(cell.execution_count is None for cell in code_cells):
        raise RuntimeError("Execute every code cell before exporting the notebook.")
    image_count = 0
    for cell in code_cells:
        for output in cell.outputs:
            if output.output_type == "error":
                raise RuntimeError(f"Notebook contains an error: {output.ename}: {output.evalue}")
            if "image/png" in output.get("data", {}):
                with Image.open(io.BytesIO(base64.b64decode(output.data["image/png"]))) as image:
                    image.verify()
                image_count += 1
    exporter = HTMLExporter(template_name="classic", embed_images=True)
    exporter.mathjax_url = ""
    exporter.require_js_url = ""
    exporter.jquery_url = ""
    exporter.raw_template = (
        '{% extends "classic/index.html.j2" %}\n'
        '{% block html_head_js %}{% endblock html_head_js %}\n'
        '{% block notebook_css %}{{ super() }}<style>{{ resources.notebook_css }}</style>{% endblock notebook_css %}\n'
    )
    resources = {"metadata": {"name": "Colorectal adenoma-carcinoma progression | Hiu Mai",
                               "path": str(root)},
                 "notebook_css": (root / "report" / "notebook.css").read_text()}
    html, _ = exporter.from_notebook_node(notebook, resources=resources)
    docs = root / "docs"
    docs.mkdir(exist_ok=True)
    if not args.export_only:
        nbformat.write(notebook, notebook_path)
    (docs / "index.html").write_text(html)
    (docs / ".nojekyll").touch()
    print(f"Executed code cells: {len(code_cells)}")
    print(f"Verified PNG outputs: {image_count}")
    print(f"Notebook: {notebook_path}")
    print(f"HTML: {docs / 'index.html'}")


if __name__ == "__main__":
    main()
