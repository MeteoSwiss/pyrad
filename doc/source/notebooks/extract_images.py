import json
import base64
from pathlib import Path

notebook_path = "plot_radar_data.ipynb"
output_dir = Path("images")
output_dir.mkdir(exist_ok=True)

with open(notebook_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if "attachments" in cell:
        for name, data in cell["attachments"].items():
            if "image/png" in data:
                img_data = base64.b64decode(data["image/png"])
                out_file = output_dir / name
                with open(out_file, "wb") as f:
                    f.write(img_data)
                print(f"Saved: {out_file}")
