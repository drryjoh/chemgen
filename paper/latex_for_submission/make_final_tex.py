#!python3
import re
import shutil
from pathlib import Path

# Settings
input_tex = Path("initial_communication_v2.tex")
output_tex = Path("final.tex")
output_dir = Path("renamed_images")
output_dir.mkdir(exist_ok=True)

# Read .tex file
with input_tex.open("r") as f:
    content = f.read()

# Match \includegraphics[...]{...} or \includegraphics{...}
pattern = r"(\\includegraphics(?:\[[^\]]*\])?)\{([^}]+)\}"
matches = re.findall(pattern, content)

# Process each match
for i, (prefix, img_path) in enumerate(matches):
    original = Path("..") /Path(img_path)
    if original.suffix == "":
        original = original.with_suffix(".png")
    ext = original.suffix if original.suffix else ".png"
    new_name = f"image_{i}{ext}"
    new_path = output_dir / new_name

    # Copy file if it exists
    full_path = input_tex.parent / original
    if full_path.exists():
        shutil.copy(full_path, new_path)
    else:
        print(f"Warning: {full_path} not found. Skipping copy.")

    # Update LaTeX content
    content = content.replace(f"{prefix}{{{img_path}}}", f"{prefix}{{{new_path.as_posix()}}}")

# Write final .tex file
with output_tex.open("w") as f:
    f.write(content)

print(f"Updated file written to {output_tex}")

