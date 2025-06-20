#!python3
import pyvista as pv
pv.start_xvfb()
import matplotlib.pyplot as plt
import numpy as np

vtk_file = "test_RR_71897.vtk"
fields = {
    "RROFO": {"cmap": "inferno", "label": "$R_{of}$ [Kg/m$^3$/s]", "vmin": None, "vmax": None},
    "RRCGO":  {"cmap": "inferno",  "label": "$R_{cg}$ [Kg/m$^3$/s]", "vmin": None, "vmax": None},
}


dpi = 300
scale = 1

# Load full mesh
mesh = pv.read(vtk_file)

print("[Field Data Keys]", mesh.field_data.keys())
for key in mesh.field_data.keys():
    print(f"{key}:", mesh.field_data[key])

bounds = mesh.bounds
x_min, x_max = bounds[0], bounds[1]
y_min, y_max = bounds[2], bounds[3]
z_min, z_max = bounds[4], bounds[5]

domain_width = x_max - x_min
domain_height = y_max - y_min

x_crop_min = x_min + 0.675 * domain_width
x_crop_max = x_min + 0.75 * domain_width

for field, cfg in fields.items():
    if field not in mesh.array_names:
        print(f"[Warning] Field '{field}' not found. Skipping.")
        continue

    data = mesh[field]
    vmin = cfg["vmin"] if cfg["vmin"] is not None else data.min()
    vmax = cfg["vmax"] if cfg["vmax"] is not None else data.max()
    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-8

    # Full domain render as usual
    window_height_px = 800
    window_width_px = int(window_height_px * domain_width / domain_height)

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
    p.camera.SetParallelScale(domain_height / 2 * 1.02)

    img = p.screenshot(return_img=True, scale=scale)

    # Plot using matplotlib, cropped via `extent`
    crop_width = x_crop_max - x_crop_min
    crop_height = y_max - y_min
    fig_width_in = 10
    fig_height_in = fig_width_in * (crop_height / crop_width)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    # Compute figure size
    crop_width = x_crop_max - x_crop_min
    crop_height = y_max - y_min
    fig_width_in = 10
    fig_height_in = fig_width_in * (crop_height / crop_width)

    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), dpi=dpi)

    # Show image with physical extents
    im = ax.imshow(img, extent=[x_min, x_max, y_min, y_max], aspect='auto')
    ax.set_xlim(x_crop_min, x_crop_max)
    ax.set_ylim(y_min, y_max)

    # Label axes
    ax.set_xlabel("x [m]", fontsize=16)
    ax.set_ylabel("y [m]", fontsize=16)
    ax.tick_params(labelsize=16)

    # Add colorbar
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm=norm, cmap=cfg["cmap"])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(cfg["label"], fontsize=20)
    cbar.ax.tick_params(labelsize=10)

    plt.savefig(f"{field}_cropped.png", dpi=dpi, bbox_inches='tight', pad_inches=0.05)
    plt.close()
