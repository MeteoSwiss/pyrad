#!/usr/bin/env python3

import argparse
import glob
import os
import re
from collections import defaultdict
from datetime import datetime
from PIL import Image, ImageDraw, ImageFont
import math

DATE_REGEX_1 = re.compile(r"(\d{4}-\d{2}-\d{2})")
DATE_REGEX_2 = re.compile(r"(\d{4}\d{2}\d{2})")


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


def create_composite(image_groups, date_order, title=None, title_size=None):
    from PIL import ImageFont

    date_keys = sorted(image_groups)
    grouped_images = [load_images(image_groups[date]) for date in date_keys]

    max_w = max(img.width for group in grouped_images for img in group)
    max_h = max(img.height for group in grouped_images for img in group)
    pad = 10
    pad_title = 20

    if date_order == "row":
        rows = len(grouped_images)
        cols = max(len(group) for group in grouped_images)
    else:
        cols = len(grouped_images)
        rows = max(len(group) for group in grouped_images)

    canvas_w = cols * (max_w + pad) + pad

    # Choose title font size
    auto_font_size = max(16, int(canvas_w * 0.02))  # 3% of width or minimum 16px
    font_size = title_size if title_size else auto_font_size

    font = ImageFont.load_default(font_size)

    title_h = font_size + pad_title if title else 0
    canvas_h = rows * (max_h + pad) + pad + title_h

    composite = Image.new("RGB", (canvas_w, canvas_h), color="white")
    draw = ImageDraw.Draw(composite)

    # Draw title centered
    if title:
        text_w = draw.textlength(title, font=font)
        draw.text(((canvas_w - text_w) // 2, pad), title, font=font, fill="black")

    y_offset = title_h

    for i, (date, images) in enumerate(zip(date_keys, grouped_images)):
        for j, img in enumerate(images):
            if date_order == "row":
                x = pad + j * (max_w + pad)
                y = y_offset + pad + i * (max_h + pad)
            else:
                x = pad + i * (max_w + pad)
                y = y_offset + pad + j * (max_h + pad)

            composite.paste(img.resize((max_w, max_h)), (x, y))

    return composite


def main():
    args = parse_args()
    image_groups = collect_images(args.patterns)

    if not image_groups:
        print("No images found.")
        return

    composite = create_composite(
        image_groups, args.date_order, title=args.title, title_size=args.title_size
    )
    composite.save(args.output)
    print(f"Composite image saved to {args.output}")


if __name__ == "__main__":
    main()
