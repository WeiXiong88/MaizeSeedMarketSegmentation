from pathlib import Path
import os
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont



OBS_FILE = "../Data/For Wei_MaizeTrials_in_Africa.xlsx"
MODEL_FILES = {
    "EPIC": "../Results/MaizeTrail_Maturity_Sim_EPIC.csv",
    "CERES-Maize": "../Results/MaizeTrail_Maturity_Sim_MZCER.csv",
    "CSM-IXIM": "../Results/MaizeTrail_Maturity_Sim_MZIXM.csv",
}

CLASS_ORDER = ["EHYB", "IHYB", "LHYB"]
SOURCE_ORDER = ["Observed", "EPIC", "CERES-Maize", "CSM-IXIM"]
MEGAENV_ORDER = [
    "Wet Upper Mid-altitude",
    "Dry Lowland",
    "Dry Mid-altitude",
    "Wet Lower Mid-altitude",
    "Highland",
    "Wet Lowland",
]

SOURCE_COLORS = {
    "Observed": "#D0D0D0",
    "EPIC": "#4E79A7",
    "CERES-Maize": "#F28E2B",
    "CSM-IXIM": "#59A14F",
}
YEAR_COLORS = {
    2011: "#4E79A7",
    2012: "#59A14F",
    2013: "#F28E2B",
    2014: "#E15759",
    2015: "#B07AA1",
    2016: "#8CD17D",
}


def hex_to_rgb(color):
    color = color.lstrip("#")
    return tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))


def load_font(size, bold=False):
    #Use regular Arial for all figure text to keep font wright consistent
    font_names = ["arial.rrf", "calibri.ttf"]
    for name in font_names:
        font_path = Path("C:/Windows/Fonts") / name
        if font_path.exists():
            return ImageFont.truetype(str(font_path), size=size)
    return ImageFont.load_default()


def normalize_megaenv(value):
    # Harmonize a small spelling variant in the trial workbook.
    if pd.isna(value):
        return np.nan
    value = str(value).strip()
    replacements = {
        "Wet Upper Mid-altitud": "Wet Upper Mid-altitude",
    }
    return replacements.get(value, value)


def load_observed_yield():
    # Oberved yield uses GWT when present and falls back to FWT
    obs = pd.read_excel(
        OBS_FILE,
        sheet_name="in",
        usecols=["Year", "Class", "MegaEnviron", "GWT", "FWT"],
    )
    obs["Yield"] = obs["GWT"].where(obs["GWT"].notna(), obs["FWT"])
    obs = obs.loc[obs["Class"].isin(CLASS_ORDER), ["Year", "Class", "MegaEnviron", "Yield"]].copy()
    obs["Yield"] = pd.to_numeric(obs["Yield"], errors="coerce")
    obs["Year"] = pd.to_numeric(obs["Year"], errors="coerce").astype("Int64")
    obs["MegaEnviron"] = obs["MegaEnviron"].map(normalize_megaenv)
    obs = obs.replace([np.inf, -np.inf], np.nan).dropna(subset=["Yield", "Year"])
    obs = obs.loc[(obs["Yield"] >= 0) & (obs["Yield"] <= 17)].copy()
    obs["Source"] = "Observed"
    return obs


def load_model_yield():
    # Stack the three crop-model outputs into one long table for plotting.
    records = []
    for model_name, file_path in MODEL_FILES.items():
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Missing model output: {file_path}")
        df = pd.read_csv(file_path)
        missing = {"Class", "MegaEnviron", "sim_yld"} - set(df.columns)
        if missing:
            raise ValueError(f"{file_path.name} is missing columns: {sorted(missing)}")
        df = df.loc[df["Class"].isin(CLASS_ORDER), ["Class", "MegaEnviron", "sim_yld"]].copy()
        df["Yield"] = pd.to_numeric(df["sim_yld"], errors="coerce")
        df["MegaEnviron"] = df["MegaEnviron"].map(normalize_megaenv)
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["Yield"])
        df = df.loc[df["Yield"] >= 0, ["Class", "MegaEnviron", "Yield"]].copy()
        df["Source"] = model_name
        records.append(df)
    return pd.concat(records, ignore_index=True)


def write_summary(observed, simulated):
    combined = pd.concat([observed, simulated], ignore_index=True)
    combined["Panel"] = combined["MegaEnviron"].fillna("All")
    all_panel = combined.copy()
    all_panel["Panel"] = "All"
    summary_input = pd.concat([all_panel, combined.dropna(subset=["MegaEnviron"])], ignore_index=True)
    summary = (
        summary_input.groupby(["Panel", "Class", "Source"], observed=True)["Yield"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            q05=lambda x: x.quantile(0.05),
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
            q95=lambda x: x.quantile(0.95),
            minimum="min",
            maximum="max",
        )
        .reset_index()
    )
    summary.to_csv("../Plots/Supp_S5_yield_maturity_model_summary.csv", index=False)


def box_stats(values):
    # Whiskers use the 5th and 95th percentiles, matching the fiture note
    return dict(zip(["q05", "q25", "q50", "q75", "q95"], np.quantile(values, [0.05, 0.25, 0.5, 0.75, 0.95])))


def draw_centered_text(draw, xy, text, font, fill=(0, 0, 0), anchor="mm"):
    draw.text(xy, text, font=font, fill=fill, anchor=anchor)


def draw_rotated_label(base, xy, text, font, fill=(0, 0, 0)):
    # PIL rotates text by drawing it on a tranparent laryer first.
    bbox = font.getbbox(text)
    width = bbox[2] - bbox[0] + 30
    height = bbox[3] - bbox[1] + 30
    label = Image.new("RGBA", (width, height), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label)
    label_draw.text((15, 15), text, font=font, fill=fill)
    rotated = label.rotate(90, expand=True)
    base.alpha_composite(rotated, (int(xy[0] - rotated.width / 2), int(xy[1] - rotated.height / 2)))


class PanelScale:
    # Convert panel-relative x position and yield values to image pixels
    def __init__(self, left, top, right, bottom, y_max):
        self.left = left
        self.top = top
        self.right = right
        self.bottom = bottom
        self.y_max = y_max

    def x(self, value):
        return int(self.left + value * (self.right - self.left))

    def y(self, value):
        return int(self.bottom - (value / self.y_max) * (self.bottom - self.top))


def draw_box(draw, scale, values, x_center, box_width, color, fill_box=True):
    # Draw one custom boxplot at a specified x position
    if len(values) == 0:
        return
    stats = box_stats(values)
    cx = scale.x(x_center)
    half_w = int((scale.x(x_center + box_width / 2) - scale.x(x_center - box_width / 2)) / 2)
    q05, q25, q50, q75, q95 = [scale.y(stats[k]) for k in ["q05", "q25", "q50", "q75", "q95"]]
    edge = (35, 35, 35)
    fill = hex_to_rgb(color) if fill_box else None
    draw.line((cx, q95, cx, q75), fill=edge, width=1)
    draw.line((cx, q25, cx, q05), fill=edge, width=1)
    draw.line((cx - half_w // 2, q95, cx + half_w // 2, q95), fill=edge, width=1)
    draw.line((cx - half_w // 2, q05, cx + half_w // 2, q05), fill=edge, width=1)
    draw.rectangle((cx - half_w, q75, cx + half_w, q25), fill=fill, outline=edge, width=1)
    draw.line((cx - half_w, q50, cx + half_w, q50), fill=(10, 10, 10), width=1)


def draw_panel(base, draw, rect, panel_label, panel_title, observed, simulated, y_max, show_y_label=False):
    # Draw one subplot: Observed genotype-year-site points plust observed/model boxplots
    left, top, right, bottom = rect
    scale = PanelScale(left + 70, top + 55, right - 22, bottom - 70, y_max)
    font_panel = load_font(55)
    font_title = load_font(55)
    font_tick = load_font(50)
    font_class = load_font(50)

    draw.text((left, top - 20), panel_label, font=font_panel, fill=(0, 0, 0))
    draw.text((left + 42, top - 18), panel_title, font=font_title, fill=(0, 0, 0))

    tick_step = 4 if y_max > 16 else 2
    tick_len = 12
    for tick in np.arange(0, y_max + 0.1, tick_step):
        y = scale.y(tick)
        draw.line((scale.left - tick_len, y, scale.left, y), fill=(35, 35, 35), width=3)
        #draw.line((scale.left, y, scale.right, y), fill=(222, 222, 222), width=1)
        draw_centered_text(draw, (scale.left - 30, y), f"{tick:g}", font_tick, anchor="rm")
    #draw.line((scale.left, scale.bottom, scale.right, scale.bottom), fill=(35, 35, 35), width=3)
    #draw.line((scale.left, scale.top, scale.left, scale.bottom), fill=(35, 35, 35), width=3)
    draw.rectangle((scale.left, scale.top, scale.right, scale.bottom), outline=(35, 35, 35), width=3)

    class_centers = np.array([0.19, 0.50, 0.81])
    offsets = {"Observed": -0.065, "EPIC": 0.015, "CERES-Maize": 0.058, "CSM-IXIM": 0.101}
    box_widths = {"Observed": 0.065, "EPIC": 0.030, "CERES-Maize": 0.030, "CSM-IXIM": 0.030}
    for center in class_centers:
        x = scale.x(center + 0.02)
        draw.line((x, scale.bottom, x, scale.bottom + tick_len), fill=(35, 35, 35), width=3)
    jitter_seed = sum((i + 1) * ord(char) for i, char in enumerate(panel_title))
    rng = np.random.default_rng(jitter_seed)

    #scatter_layer = Image.new("RGBA", base.size, (255, 255, 255, 0))
    #scatter_draw = ImageDraw.Draw(scatter_layer)
    #for center, maturity_class in zip(class_centers, CLASS_ORDER):
    #    obs_class = observed.loc[observed["Class"] == maturity_class]
    #    x_center = center + offsets["Observed"]
    #    jitter = rng.normal(0, 0.011, size=len(obs_class))
    #    for x_val, (_, row) in zip(x_center + jitter, obs_class.iterrows()):
    #        year = int(row["Year"])
    #        color = hex_to_rgb(YEAR_COLORS.get(year, "#777777"))
    #        px, py = scale.x(x_val), scale.y(row["Yield"])
    #        scatter_draw.ellipse((px - 1, py - 1, px + 1, py + 1), fill=(*color, 38))  #Changed 3 to 1
    #base.alpha_composite(scatter_layer)

    for center, maturity_class in zip(class_centers, CLASS_ORDER):
        obs_values = observed.loc[observed["Class"] == maturity_class, "Yield"].to_numpy()
        #draw_box(draw, scale, obs_values, center + offsets["Observed"], box_widths["Observed"], SOURCE_COLORS["Observed"])
        draw_box(
            draw,
            scale,
            obs_values,
            center + offsets["Observed"],
            box_widths["Observed"],
            SOURCE_COLORS["Observed"],
            fill_box=False,
        )
        for source in ["EPIC", "CERES-Maize", "CSM-IXIM"]:
            sim_values = simulated.loc[
                (simulated["Class"] == maturity_class) & (simulated["Source"] == source),
                "Yield",
            ].to_numpy()
            draw_box(draw, scale, sim_values, center + offsets[source], box_widths[source], SOURCE_COLORS[source])
        draw_centered_text(draw, (scale.x(center + 0.02), scale.bottom + 45), maturity_class, font_class)

    scatter_layer = Image.new("RGBA", base.size, (255, 255, 255, 0))
    scatter_draw = ImageDraw.Draw(scatter_layer)
    for center, maturity_class in zip(class_centers, CLASS_ORDER):
        obs_class = observed.loc[observed["Class"] == maturity_class]
        x_center = center + offsets["Observed"]
        jitter = rng.normal(0, 0.011, size=len(obs_class))
        for x_val, (_, row) in zip(x_center + jitter, obs_class.iterrows()):
            year = int(row["Year"])
            color = hex_to_rgb(YEAR_COLORS.get(year, "#777777"))
            px, py = scale.x(x_val), scale.y(row["Yield"])
            scatter_draw.ellipse((px - 3, py - 3, px + 3, py + 3), fill=(*color, 38))  #Changed 3 to 1
    base.alpha_composite(scatter_layer)

    if show_y_label:
        draw_rotated_label(base, (left - 45, (scale.top + scale.bottom) / 2), "Maize yield (t/ha)", load_font(55))
        #                                40                                                                    26  


def draw_legend_panel(draw, rect, years):
    # The empty eight panel is reserved for source and trial-year legends
    left, top, right, bottom = rect
    font_panel = load_font(55)
    font_heading = load_font(50)
    font_text = load_font(50)
    draw.text((left, top + 4), "h", font=font_panel, fill=(0, 0, 0))

    x0, y0 = left + 55, top + 65
    draw.text((x0, y0), "Yield source", font=font_heading, fill=(0, 0, 0))
    y = y0 + 50
    for source in SOURCE_ORDER:
        fill = None if source == "Observed" else hex_to_rgb(SOURCE_COLORS[source])
        draw.rectangle((x0, y, x0 + 45, y + 28), fill=fill, outline=(35, 35, 35), width=2)
        draw.text((x0 + 62, y - 4), source, font=font_text, fill=(0, 0, 0))
        y += 43

    y += 24
    draw.text((x0, y), "Observed trial year", font=font_heading, fill=(0, 0, 0))
    y += 50
    for year in years:
        color = hex_to_rgb(YEAR_COLORS.get(int(year), "#777777"))
        draw.ellipse((x0 + 8, y + 3, x0 + 34, y + 29), fill=color, outline=(35, 35, 35), width=1)
        draw.text((x0 + 62, y - 3), str(int(year)), font=font_text, fill=(0, 0, 0))
        y += 42

    y += 22
    note = [
        "Wide boxes: observed trials",
        "Narrow boxes: crop models",
        "Whiskers: 5th-95th percentile",
    ]
    for line in note:
        draw.text((x0, y), line, font=font_text, fill=(35, 35, 35))
        y += 38


def plot_s5():
    observed = load_observed_yield()
    simulated = load_model_yield()
    write_summary(observed, simulated)

    width, height = 3600, 2050
    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image)

    all_values = pd.concat([observed["Yield"], simulated["Yield"]], ignore_index=True)
    y_max = 17 #float(np.ceil(all_values.max() / 2.0) * 2.0)
    #y_max = max(y_max, 12.0)

    margin_left, margin_top = 95, 65 #115, 85
    panel_w, panel_h = 820, 850 #810, 840
    gap_x, gap_y = 45, 95 #55, 115
    panels = []
    for row in range(2):
        for col in range(4):
            left = margin_left + col * (panel_w + gap_x)
            top = margin_top + row * (panel_h + gap_y)
            panels.append((left, top, left + panel_w, top + panel_h))

    panel_specs = [("a", "All trial sites", None)] + [
        (chr(ord("b") + i), megaenv, megaenv) for i, megaenv in enumerate(MEGAENV_ORDER)
    ]
    for idx, (label, title, megaenv) in enumerate(panel_specs):
        if megaenv is None:
            obs_panel = observed
            sim_panel = simulated
        else:
            obs_panel = observed.loc[observed["MegaEnviron"] == megaenv]
            sim_panel = simulated.loc[simulated["MegaEnviron"] == megaenv]
        draw_panel(
            image,
            draw,
            panels[idx],
            label,
            title,
            obs_panel,
            sim_panel,
            y_max,
            show_y_label=idx in [0, 4],
        )

    draw_legend_panel(draw, panels[-1], sorted(observed["Year"].dropna().unique()))
    draw_centered_text(
        draw,
        (width / 2 - 180, height - 90),
        "Maturity class",
        load_font(60),
    )

    output_png = "../Plots/plot_s5.png"
    output_pdf = "../Plots/plot_s5.pdf"
    image_rgb = image.convert("RGB")
    image_rgb.resize((1800, 1025), Image.Resampling.LANCZOS).save(output_png, dpi=(800, 800))
    image_rgb.save(output_pdf, "PDF", resolution=800)


if __name__ == "__main__":
    plot_s5()
