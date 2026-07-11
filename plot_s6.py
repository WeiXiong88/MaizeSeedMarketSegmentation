from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "Data"
RESULTS_DIR = SCRIPT_DIR.parent / "Results"
PLOTS_DIR = SCRIPT_DIR.parent / "Plots"

BENCHMARK_FILES = {
    "EPIC": RESULTS_DIR / "00_simunit_05degree_EPIC__regionalSim_withBenchmarkCultivar.csv",
    "CERES-Maize": RESULTS_DIR / "00_simunit_05degree_MZCER_regionalSim_withBenchmarkCultivar.csv",
    "CSM-IXIM": RESULTS_DIR / "00_simunit_05degree_MZIXM_regionalSim_withBenchmarkCultivar.csv",
}
CALI_FILES = {
    "EPIC": RESULTS_DIR / "00_simunit_05degree_EPIC__withCaliPara.csv",
    "CERES-Maize": RESULTS_DIR / "00_simunit_05degree_MZCER__withCaliPara.csv",
    "CSM-IXIM": RESULTS_DIR / "00_simunit_05degree_MZIXM__withCaliPara.csv",
}
SIM5_FILE = DATA_DIR / "01_simunit_5min_EPIC_4_Simulation.csv"
HVSTAT_FILE = DATA_DIR / "hvstat_africa_data_v1.0.csv"

OUTPUT_PNG = PLOTS_DIR / "plot_s6.png"
OUTPUT_PDF = PLOTS_DIR / "plot_s6.pdf"
STATS_CSV = RESULTS_DIR / "plot_s6_panel_statistics.csv"
PROD_MATCH_CSV = RESULTS_DIR / "plot_s6_subnational_production_matches.csv"

MODEL_ORDER = ["EPIC", "CERES-Maize", "CSM-IXIM"]
YEAR_COLUMNS = [f"sy{year:02d}" for year in range(6, 16)]


def load_font(size):
    """Use regular Arial when available, with Calibri as a Windows fallback."""
    for name in ["arial.ttf", "calibri.ttf"]:
        path = Path("C:/Windows/Fonts") / name
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def normalize_name(value):
    """Normalize administrative names for robust joins across grid and harvest-stat files."""
    if pd.isna(value):
        return "none"
    text = str(value).strip().lower()
    return "none" if text in {"", "nan", "none", "na", "n/a"} else text


def clean_xy(df, x_col, y_col):
    """Keep finite, positive observed and simulated values only."""
    out = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    out = out[(out[x_col] > 0) & (out[y_col] > 0)]
    return out.rename(columns={x_col: "obs", y_col: "sim"})


def load_benchmark_grid(model):
    """Row 1: compare SPAM/old observed yield obsy1 with benchmark-cultivar simulations."""
    df = pd.read_csv(BENCHMARK_FILES[model])
    df["sim_mean_2006_2015"] = df[YEAR_COLUMNS].replace(0, np.nan).mean(axis=1, skipna=True)
    return clean_xy(df, "obsy1", "sim_mean_2006_2015")


def load_calibrated_grid(model):
    """Row 2: compare yld_se1 observation with calibrated simulated yield."""
    df = pd.read_csv(CALI_FILES[model], usecols=["lon", "lat", "yld_se1", "cali_yld_se1"])
    return clean_xy(df, "yld_se1", "cali_yld_se1")


def load_hvstat_production():
    """Mean recorded maize production for 2005-2015 from the harvest-stat data."""
    hv = pd.read_csv(HVSTAT_FILE)
    hv = hv[
        (hv["product"].str.lower() == "maize")
        & (hv["harvest_year"].between(2005, 2015))
    ].copy()
    hv["admin_2"] = hv["admin_2"].fillna("none")
    hv["level"] = np.where(hv["admin_2"].map(normalize_name).eq("none"), "adm1", "adm2")
    hv = (
        hv.groupby(["country", "admin_1", "admin_2", "level"], dropna=False)["production"]
        .mean()
        .reset_index()
    )
    hv["k0"] = hv["country"].map(normalize_name)
    hv["k1"] = hv["admin_1"].map(normalize_name)
    hv["k2"] = hv["admin_2"].map(normalize_name)
    return hv


def load_calibrated_production(model, hv_obs):
    """Row 3: aggregate calibrated grid yield to subnational production and match harvest-stat records."""
    sim5 = pd.read_csv(SIM5_FILE, usecols=["lon", "lat", "adm0", "adm1", "adm2", "mzarea"])
    yld = pd.read_csv(CALI_FILES[model], usecols=["lon", "lat", "cali_yld_se1"])
    df = sim5.merge(yld, on=["lon", "lat"], how="left")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["cali_yld_se1", "mzarea"])
    df = df[(df["cali_yld_se1"] > 0) & (df["mzarea"] > 0)].copy()
    df["sim_production"] = df["cali_yld_se1"] * df["mzarea"]
    df["k0"] = df["adm0"].map(normalize_name)
    df["k1"] = df["adm1"].map(normalize_name)
    df["k2"] = df["adm2"].map(normalize_name)

    # Harvest-stat admin_2 == none is an admin-1 total, so compare it with admin-1 grid sums.
    sim_adm1 = df.groupby(["k0", "k1"], dropna=False)["sim_production"].sum().reset_index()
    obs_adm1 = hv_obs[hv_obs["level"] == "adm1"].copy()
    matched_adm1 = obs_adm1.merge(sim_adm1, on=["k0", "k1"], how="inner")

    # Harvest-stat rows with a named admin_2 are compared with admin-2 grid sums.
    sim_adm2 = df.groupby(["k0", "k1", "k2"], dropna=False)["sim_production"].sum().reset_index()
    obs_adm2 = hv_obs[hv_obs["level"] == "adm2"].copy()
    matched_adm2 = obs_adm2.merge(sim_adm2, on=["k0", "k1", "k2"], how="inner")

    matched = pd.concat([matched_adm1, matched_adm2], ignore_index=True, sort=False)
    matched["model"] = model
    matched["obs"] = matched["production"] / 1_000_000.0
    matched["sim"] = matched["sim_production"] / 1_000_000.0
    return clean_xy(matched, "obs", "sim"), matched


def panel_stats(obs, sim):
    """Compute compact validation metrics for each panel."""
    obs = np.asarray(obs, dtype=float)
    sim = np.asarray(sim, dtype=float)
    n = len(obs)
    if n == 0:
        return {"n": 0, "r2": np.nan, "rmse": np.nan, "mae": np.nan, "bias": np.nan}
    r = np.corrcoef(obs, sim)[0, 1] if n > 1 else np.nan
    return {
        "n": int(n),
        "r2": float(r**2) if np.isfinite(r) else np.nan,
        "rmse": float(np.sqrt(np.mean((sim - obs) ** 2))),
        "mae": float(np.mean(np.abs(sim - obs))),
        "bias": float(np.mean(sim - obs)),
    }


def color_from_density(x, y, xlim, ylim):
    """Approximate point density by assigning each point the count of its 2-D bin."""
    bins = 60
    xi = np.clip(((x - xlim[0]) / (xlim[1] - xlim[0]) * bins).astype(int), 0, bins - 1)
    yi = np.clip(((y - ylim[0]) / (ylim[1] - ylim[0]) * bins).astype(int), 0, bins - 1)
    counts = np.zeros((bins, bins), dtype=int)
    np.add.at(counts, (xi, yi), 1)
    density = counts[xi, yi].astype(float)
    if density.max() > density.min():
        density = (density - density.min()) / (density.max() - density.min())
    else:
        density[:] = 0.0
    return density


def viridis_rgb(value):
    """Small viridis-like palette interpolator for density-colored scatter points."""
    palette = np.array(
        [
            (68, 1, 84),
            (59, 82, 139),
            (33, 145, 140),
            (94, 201, 98),
            (253, 231, 37),
        ],
        dtype=float,
    )
    value = float(np.clip(value, 0, 1))
    pos = value * (len(palette) - 1)
    lo = int(np.floor(pos))
    hi = min(lo + 1, len(palette) - 1)
    frac = pos - lo
    return tuple(np.round(palette[lo] * (1 - frac) + palette[hi] * frac).astype(int))


class Panel:
    """Map data coordinates to one subplot's pixel coordinates."""

    def __init__(self, left, top, right, bottom, xlim, ylim):
        self.left = left
        self.top = top
        self.right = right
        self.bottom = bottom
        self.xlim = xlim
        self.ylim = ylim

    def x(self, value):
        return int(self.left + (value - self.xlim[0]) / (self.xlim[1] - self.xlim[0]) * (self.right - self.left))

    def y(self, value):
        return int(self.bottom - (value - self.ylim[0]) / (self.ylim[1] - self.ylim[0]) * (self.bottom - self.top))


def draw_text_center(draw, xy, text, font, fill=(0, 0, 0), anchor="mm"):
    draw.text(xy, text, font=font, fill=fill, anchor=anchor)


def draw_rotated_label(image, xy, text, font):
    bbox = font.getbbox(text)
    label = Image.new("RGBA", (bbox[2] - bbox[0] + 40, bbox[3] - bbox[1] + 40), (255, 255, 255, 0))
    d = ImageDraw.Draw(label)
    d.text((20, 20), text, font=font, fill=(0, 0, 0))
    rotated = label.rotate(90, expand=True)
    image.alpha_composite(rotated, (int(xy[0] - rotated.width / 2), int(xy[1] - rotated.height / 2)))


def draw_panel(image, draw, rect, data, model, panel_letter, x_label, y_label, xlim, ylim):
    """Draw one observed-vs-simulated validation scatter panel."""
    left, top, right, bottom = rect
    panel = Panel(left + 92, top + 70, right - 35, bottom - 92, xlim, ylim)
    font_letter = load_font(44)
    font_title = load_font(34)
    font_tick = load_font(26)
    font_axis = load_font(30)
    font_stats = load_font(24)

    x = data["obs"].to_numpy(dtype=float)
    y = data["sim"].to_numpy(dtype=float)
    keep = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    x_plot = x[keep]
    y_plot = y[keep]

    draw.text((left, top + 8), panel_letter, font=font_letter, fill=(0, 0, 0))
    draw.text((left + 60, top + 18), model, font=font_title, fill=(0, 0, 0))
    draw.rectangle((panel.left, panel.top, panel.right, panel.bottom), outline=(0, 0, 0), width=3)

    tick_step = 2 if xlim[1] <= 14 else 1
    for tick in np.arange(xlim[0], xlim[1] + 0.001, tick_step):
        px = panel.x(tick)
        draw.line((px, panel.bottom, px, panel.bottom + 12), fill=(0, 0, 0), width=3)
        draw_text_center(draw, (px, panel.bottom + 43), f"{tick:g}", font_tick)
    for tick in np.arange(ylim[0], ylim[1] + 0.001, tick_step):
        py = panel.y(tick)
        draw.line((panel.left - 12, py, panel.left, py), fill=(0, 0, 0), width=3)
        draw_text_center(draw, (panel.left - 24, py), f"{tick:g}", font_tick, anchor="rm")

    # 1:1 line
    lo = max(xlim[0], ylim[0])
    hi = min(xlim[1], ylim[1])
    draw.line((panel.x(lo), panel.y(lo), panel.x(hi), panel.y(hi)), fill=(220, 35, 35), width=4)

    densities = color_from_density(x_plot, y_plot, xlim, ylim) if len(x_plot) else []
    order = np.argsort(densities) if len(x_plot) else []
    for idx in order:
        px = panel.x(x_plot[idx])
        py = panel.y(y_plot[idx])
        color = viridis_rgb(densities[idx])
        draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=(*color, 175))

    stats = panel_stats(x, y)
    stats_text = (
        f"n = {stats['n']}\n"
        f"R2 = {stats['r2']:.2f}\n"
        f"RMSE = {stats['rmse']:.2f}\n"
        f"Bias = {stats['bias']:.2f}"
    )
    box_left, box_top = panel.left + 18, panel.top + 18
    draw.rounded_rectangle(
        (box_left - 8, box_top - 8, box_left + 165, box_top + 132),
        radius=6,
        fill=(255, 255, 255, 225),
        outline=(210, 210, 210),
        width=2,
    )
    draw.multiline_text((box_left, box_top), stats_text, font=font_stats, fill=(0, 0, 0), spacing=7)

    draw_text_center(draw, ((panel.left + panel.right) / 2, bottom - 22), x_label, font_axis)
    draw_rotated_label(image, (left + 24, (panel.top + panel.bottom) / 2), y_label, font_axis)
    return stats


def build_plot():
    """Build the new Supplementary Fig. S6 validation figure and supporting statistics."""
    rows = []
    plot_data = {}

    for model in MODEL_ORDER:
        data = load_benchmark_grid(model)
        plot_data[("benchmark", model)] = data
        rows.append({"row": "benchmark cultivar grid yield", "model": model, **panel_stats(data["obs"], data["sim"])})

    for model in MODEL_ORDER:
        data = load_calibrated_grid(model)
        plot_data[("calibrated", model)] = data
        rows.append({"row": "calibrated grid yield", "model": model, **panel_stats(data["obs"], data["sim"])})

    hv_obs = load_hvstat_production()
    matched_outputs = []
    for model in MODEL_ORDER:
        data, matched = load_calibrated_production(model, hv_obs)
        plot_data[("production", model)] = data
        matched_outputs.append(matched)
        rows.append({"row": "subnational HVSTAT production", "model": model, **panel_stats(data["obs"], data["sim"])})

    pd.DataFrame(rows).to_csv(STATS_CSV, index=False)
    pd.concat(matched_outputs, ignore_index=True).to_csv(PROD_MATCH_CSV, index=False)

    width, height = 3000, 3000
    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image)
    font_row = load_font(32)

    margin_left, margin_top = 120, 95
    panel_w, panel_h = 900, 860
    gap_x, gap_y = 40, 82
    letters = list("abcdefghi")
    row_specs = [
        ("benchmark", "Benchmark cultivar grid yield", "Observed yield, SPAM-obsy1 (t ha-1)", "Simulated yield, 2006-2015 mean (t ha-1)", (0, 14), (0, 14)),
        ("calibrated", "Grid-calibrated cultivar yield", "Observed yield, yld_se1 (t ha-1)", "Simulated yield, calibrated (t ha-1)", (0, 14), (0, 14)),
        ("production", "Subnational production validation", "Observed production, 2005-2015 mean (million t)", "Simulated production (million t)", (0, 4), (0, 4)),
    ]

    letter_i = 0
    for r, (row_key, row_title, x_label, y_label, xlim, ylim) in enumerate(row_specs):
        row_top = margin_top + r * (panel_h + gap_y)
        draw.text((margin_left, row_top - 58), row_title, font=font_row, fill=(0, 0, 0))
        for c, model in enumerate(MODEL_ORDER):
            left = margin_left + c * (panel_w + gap_x)
            top = row_top
            rect = (left, top, left + panel_w, top + panel_h)
            draw_panel(
                image,
                draw,
                rect,
                plot_data[(row_key, model)],
                model,
                letters[letter_i],
                x_label,
                y_label,
                xlim,
                ylim,
            )
            letter_i += 1

    image_rgb = image.convert("RGB")
    image_rgb.save(OUTPUT_PNG, dpi=(300, 300))
    image_rgb.save(OUTPUT_PDF, "PDF", resolution=300)
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_PDF}")
    print(f"Wrote {STATS_CSV}")
    print(f"Wrote {PROD_MATCH_CSV}")


if __name__ == "__main__":
    build_plot()
