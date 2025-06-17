#!python3
import pyvista as pv
pv.start_xvfb()
import matplotlib.pyplot as plt
import numpy as np

# File
vtk_file = "VTK/yass_with_fouier_ICs_71246.vtk"

# Field configs
full_field = {
    "maxp": {"cmap": "coolwarm", "label": "Maximum Pressure", "vmin": None, "vmax": 300000},
}

dpi = 300
scale = 2

# Load mesh
mesh = pv.read(vtk_file)
bounds = mesh.bounds
x_min, x_max = bounds[0], bounds[1]
x_50 = x_min + 0.5 * (x_max - x_min)
x_80 = x_min + 0.8 * (x_max - x_min)
y_min, y_max = bounds[2], bounds[3]
z_min, z_max = bounds[4], bounds[5]
height = y_max - y_min
width = x_max - x_min
sub_width = x_80 - x_50

# --------------------------
# 1. Full-domain maxp image
# --------------------------
for field, cfg in full_field.items():
    if field not in mesh.array_names:
        print(f"[Warning] Field '{field}' not found. Skipping.")
        continue

    data = mesh[field]
    vmin = cfg["vmin"] if cfg["vmin"] is not None else data.min()
    vmax = cfg["vmax"] if cfg["vmax"] is not None else data.max()

    window_height_px = 800
    window_width_px = int(window_height_px * width / height)

    p = pv.Plotter(off_screen=True, window_size=(window_width_px, window_height_px))
    p.add_mesh(mesh, scalars=field, cmap=cfg["cmap"], clim=(vmin, vmax), show_scalar_bar=False)
    p.view_xy()
    p.camera.SetParallelProjection(True)

    xmid = 0.5 * (x_min + x_max)
    ymid = 0.5 * (y_min + y_max)
    zmid = 0.5 * (z_min + z_max)
    p.camera.SetFocalPoint(xmid, ymid, zmid)
    p.camera.SetPosition(xmid, ymid, zmid + 1)
    p.camera.SetViewUp(0, 1, 0)
    p.camera.SetParallelScale(height / 2 * 1.02)

    img = p.screenshot(return_img=True, scale=scale)
    img_height, img_width = img.shape[:2]
    figsize = (img_width / dpi, img_height / dpi)
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    # Plot with axis labels and colorbar
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Physical extents for imshow
    im = ax.imshow(img, extent=[x_min, x_max, y_min, y_max], aspect='auto')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Axis labels and ticks
    ax.set_xlabel("x [m]", fontsize=16)
    ax.set_ylabel("y [m]", fontsize=16)
    ax.tick_params(labelsize=16)

    # Colorbar setup
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm=norm, cmap=cfg["cmap"])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    from matplotlib.ticker import FuncFormatter

    # Scale colorbar tick labels
    atm_scale = 101325.0
    def format_atm(x, pos):
        return f"{x / atm_scale:.2f}"

    cbar.set_label(f"{cfg['label']} [atm]", fontsize=16)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(format_atm))


    # Save figure
    plt.savefig("maxp_full.png", dpi=dpi, bbox_inches='tight', pad_inches=0.05)
    plt.close()
