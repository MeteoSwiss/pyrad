#!/usr/bin/env python3

import argparse
import glob
import os
import re
from collections import defaultdict
from datetime import datetime
from PIL import Image, ImageDraw, ImageFont
import math

DATE_REGEX_1 = re.compile(r"(\d{4}\d{2}\d{2}\d{2}\d{2}\d{2})")
DATE_REGEX_2 = re.compile(r"(\d{4}-\d{2}-\d{2})")
DATE_REGEX_3 = re.compile(r"(\d{4}\d{2}\d{2})")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create composite image from images with date-encoded filenames."
    )
    parser.add_argument(
        "patterns", nargs="+", help="One or more glob patterns (e.g. 'figs/**/*.png')."
    )
    parser.add_argument(
        "--date-order",
        choices=["row", "column"],
        required=True,
        help="Group images with same date in a row or column.",
    )
    parser.add_argument(
        "--title", help="Optional title to place at the top of the composite image."
    )
    parser.add_argument(
        "--title-size",
        type=int,
        help="Font size for the title (in pixels). Defaults to ~3%% of image width.",
    )
    parser.add_argument(
        "--coltitles",
        help="Optional title to place at the top of every column. Separate titles should be comma-separated.",
    )
    parser.add_argument(
        "--coltitles-size",
        type=int,
        help="Font size for the column titles (in pixels). Defaults to ~2%% of image width.",
    )
    parser.add_argument(
        "--rowtitles",
        help="Optional title to place at the top of every row. Separate titles should be comma-separated.",
    )
    parser.add_argument(
        "--rowtitles-size",
        type=int,
        help="Font size for the row titles (in pixels). Defaults to ~2%% of image width.",
    )
    parser.add_argument(
        "--resize",
        type=int,
        help="Resize factor of the final image (in percent), default is 100%",
    )
    parser.add_argument(
        "--output",
        default="composite.jpg",
        help="Filename of the output composite image.",
    )

    return parser.parse_args()


def extract_date_from_filename(filename):
    match = DATE_REGEX_1.search(filename)
    if match:
        return match.group(1)
    match = DATE_REGEX_2.search(filename)
    if match:
        return match.group(1)
    match = DATE_REGEX_3.search(filename)
    if match:
        return match.group(1)
    return "unknown"


def collect_images(patterns):
    image_groups = defaultdict(list)
    for pattern in patterns:
        for path in glob.glob(pattern, recursive=True):
            if os.path.isfile(path) and path.lower().endswith(
                (".png", ".jpg", ".jpeg", ".bmp", ".gif")
            ):
                date = extract_date_from_filename(os.path.basename(path))
                image_groups[date].append(path)
    return dict(sorted(image_groups.items()))


def load_images(image_paths):
    return [Image.open(p).convert("RGB") for p in image_paths]


def create_composite(
    image_groups,
    date_order,
    title=None,
    title_size=None,
    col_titles=None,
    coltitle_size=None,
    row_titles=None,
    rowtitle_size=None,
    resize_percent=100,
):
    date_keys = sorted(image_groups)
    grouped_images = [load_images(image_groups[date]) for date in date_keys]

    max_w = max(img.width for group in grouped_images for img in group)
    max_h = max(img.height for group in grouped_images for img in group)
    pad = 10
    pad_title = 20
    # Determine layout
    if date_order == "row":
        rows = len(grouped_images)
        cols = max(len(group) for group in grouped_images)
    else:
        cols = len(grouped_images)
        rows = max(len(group) for group in grouped_images)

    # Font setup
    canvas_w = cols * (max_w + pad) + pad
    auto_font_size = max(16, int(canvas_w * 0.03))
    font_size = int(title_size) if title_size else auto_font_size
    font = ImageFont.load_default(font_size)

    auto_font_size = max(14, int(canvas_w * 0.02))
    coltitle_size = int(coltitle_size) if coltitle_size else auto_font_size
    rowtitle_size = int(rowtitle_size) if rowtitle_size else auto_font_size
    rowfont = ImageFont.load_default(rowtitle_size)
    colfont = ImageFont.load_default(coltitle_size)

    draw_dummy = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    row_title_width = (
        max(draw_dummy.textlength(title, font=rowfont) for title in row_titles) + pad
        if row_titles
        else 0
    )
    col_title_height = coltitle_size + pad_title if col_titles else 0
    main_title_height = font_size + pad_title if title else 0

    # Canvas dimensions
    canvas_w += row_title_width
    canvas_w = int(canvas_w)
    canvas_h = int(rows * (max_h + pad) + pad + main_title_height + col_title_height)

    composite = Image.new("RGB", (canvas_w, canvas_h), color="white")
    draw = ImageDraw.Draw(composite)

    # Draw main title
    if title:
        text_w = draw.textlength(title, font=font)
        draw.text(((canvas_w - text_w) // 2, pad), title, font=font, fill="black")

    y_offset = main_title_height + col_title_height

    for i, (date, images) in enumerate(zip(date_keys, grouped_images)):
        for j, img in enumerate(images):
            if date_order == "row":
                row_idx, col_idx = i, j
            else:
                row_idx, col_idx = j, i

            x = int(pad + row_title_width + col_idx * (max_w + pad))
            y = int(y_offset + row_idx * (max_h + pad))

            composite.paste(img.resize((max_w, max_h)), (x, y))

    # Draw column titles
    if col_titles:
        for j in range(cols):
            col_title = col_titles[j] if j < len(col_titles) else ""
            x = pad + row_title_width + j * (max_w + pad)
            y = main_title_height
            draw.text(
                (x + max_w // 2 - draw.textlength(col_title, font=font) // 2, y),
                col_title,
                font=colfont,
                fill="black",
            )

    # Draw row titles
    if row_titles:
        for i in range(rows):
            row_title = row_titles[i] if i < len(row_titles) else ""
            x = pad
            y = y_offset + i * (max_h + pad) + max_h // 2 - font_size // 2
            draw.text((x, y), row_title, font=rowfont, fill="black")
    # Resize if needed
    if resize_percent != 100:
        new_w = int(composite.width * resize_percent / 100)
        new_h = int(composite.height * resize_percent / 100)
        composite = composite.resize((new_w, new_h))

    return composite


def main():
    args = parse_args()
    image_groups = collect_images(args.patterns)

    if not image_groups:
        print("No images found.")
        return

    if args.rowtitles:
        args.rowtitles = args.rowtitles.split(",")
    if args.coltitles:
        args.coltitles = args.coltitles.split(",")

    composite = create_composite(
        image_groups,
        args.date_order,
        title=args.title,
        title_size=args.title_size,
        row_titles=args.rowtitles,
        rowtitle_size=args.rowtitles_size,
        col_titles=args.coltitles,
        coltitle_size=args.coltitles_size,
        resize_percent=args.resize,
    )
    composite.save(args.output)
    print(f"Composite image saved to {args.output}")


if __name__ == "__main__":
    main()
