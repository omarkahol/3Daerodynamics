import sys
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from math import radians, degrees, sin, cos
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
import matplotlib.colors as colors

# Ensure parent directory is in path to import aerodynamics package
sys.path.append(str(Path(__file__).parent))

from aerodynamics import (
    WingMesher, VortexLatticeSolver, 
    NonLinearVortexLatticeSolver, AeroFoil, VortexLineTool, setup_plot_style
)

# ANSI Escape Codes for CLI styling
RESET = "\033[0m"
BOLD = "\033[1m"
CYAN = "\033[36m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
RED = "\033[31m"
BLUE = "\033[34m"

ASCII_LOGO = fr"""{CYAN}{BOLD}
   ===============================================================================================
    ____  _____           ______ _____   ____  _____ __   __ _   _          __  __ _____ _____  _____ 
   |___ \|  __ \   /\    |  ____|  __ \ / __ \|  __ \\ \ / /| \ | |   /\   |  \/  |_   _/ ____|/ ____|
     __) | |  | | /  \   | |__  | |__) | |  | | |  | |\ V / |  \| |  /  \  | \  / | | || |    | (___  
    |__ <| |  | |/ /\ \  |  __| |  _  /| |  | | |  | |  | |  | . ` | / /\ \ | |\/| | | || |     \___ \ 
    ___) | |__| / ____ \ | |____| | \ \| |__| | |__| | | |  | |\  |/ ____ \| |  | |_| || |____ ____) |
   |____/|_____/_/    \_\|______|_|  \_\\____/|_____/   |_|  |_| \_/_/    \_\_|  |_|_____\_____|_____/ 
   ===============================================================================================
                            Vortex Lattice Method Suite
   ==============================================================================================={RESET}"""

def save_vtu(filepath, mesher, panels, solver_type, decamber_angles=None, coefficients=None, profile_drag_coeff=None):
    n_span = mesher.n_span
    n_chord = mesher.n_chord
    
    y_nodes = np.linspace(-0.5 * mesher.span, 0.5 * mesher.span, n_span)
    s_nodes = np.linspace(0.0, 1.0, n_chord)
    
    points = []
    for y in y_nodes:
        for s in s_nodes:
            chord, x_le, z_le, twist_rad, z_camber = mesher._interpolate_station_properties(y, s)
            dx = s * chord
            dz = z_camber
            x_rot = dx * np.cos(twist_rad) - dz * np.sin(twist_rad)
            z_rot = dx * np.sin(twist_rad) + dz * np.cos(twist_rad)
            points.append([x_le + x_rot, y, z_le + z_rot])
            
    points = np.array(points)
    n_points = len(points)
    
    connectivity = []
    offsets = []
    types = []
    
    offset = 0
    for i in range(n_span - 1):
        for j in range(n_chord - 1):
            p1_idx = i * n_chord + j
            p2_idx = (i + 1) * n_chord + j
            p4_idx = (i + 1) * n_chord + (j + 1)
            p3_idx = i * n_chord + (j + 1)
            
            connectivity.extend([p1_idx, p2_idx, p4_idx, p3_idx])
            offset += 4
            offsets.append(offset)
            types.append(9) # VTK_QUAD
            
    n_cells = len(types)
    
    gammas = [p.gamma for p in panels]
    induced_velocities = [p.inducedVelocity for p in panels]
    forces = [p.getForce() for p in panels]
    areas = [p.surfaceArea for p in panels]
    if decamber_angles is None:
        decamber_angles = np.zeros(n_cells)
    if profile_drag_coeff is None:
        profile_drag_coeff = np.zeros(n_cells)
        
    vtk_file = ET.Element("VTKFile", type="UnstructuredGrid", version="0.1", byte_order="LittleEndian")
    ug = ET.SubElement(vtk_file, "UnstructuredGrid")
    piece = ET.SubElement(ug, "Piece", NumberOfPoints=str(n_points), NumberOfCells=str(n_cells))
    
    # Points
    pts_el = ET.SubElement(piece, "Points")
    pts_arr = ET.SubElement(pts_el, "DataArray", type="Float32", Name="Points", NumberOfComponents="3", format="ascii")
    pts_arr.text = "\n" + " ".join(f"{pt[0]:.6f} {pt[1]:.6f} {pt[2]:.6f}" for pt in points) + "\n"
    
    # Cells
    cells_el = ET.SubElement(piece, "Cells")
    conn_arr = ET.SubElement(cells_el, "DataArray", type="Int32", Name="connectivity", format="ascii")
    conn_arr.text = "\n" + " ".join(map(str, connectivity)) + "\n"
    
    off_arr = ET.SubElement(cells_el, "DataArray", type="Int32", Name="offsets", format="ascii")
    off_arr.text = "\n" + " ".join(map(str, offsets)) + "\n"
    
    types_arr = ET.SubElement(cells_el, "DataArray", type="UInt8", Name="types", format="ascii")
    types_arr.text = "\n" + " ".join(map(str, types)) + "\n"
    
    # CellData
    cd_el = ET.SubElement(piece, "CellData")
    gamma_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="gamma", format="ascii")
    gamma_arr.text = "\n" + " ".join(f"{g:.6f}" for g in gammas) + "\n"
    
    vel_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="induced_velocity", NumberOfComponents="3", format="ascii")
    vel_arr.text = "\n" + " ".join(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}" for v in induced_velocities) + "\n"
    
    force_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="force", NumberOfComponents="3", format="ascii")
    force_arr.text = "\n" + " ".join(f"{f[0]:.6f} {f[1]:.6f} {f[2]:.6f}" for f in forces) + "\n"
    
    area_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="area", format="ascii")
    area_arr.text = "\n" + " ".join(f"{a:.6f}" for a in areas) + "\n"
    
    decamber_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="decamber_angle", format="ascii")
    decamber_arr.text = "\n" + " ".join(f"{d:.6f}" for d in decamber_angles) + "\n"
    
    cdp_arr = ET.SubElement(cd_el, "DataArray", type="Float32", Name="profile_drag_coeff", format="ascii")
    cdp_arr.text = "\n" + " ".join(f"{c:.6f}" for c in profile_drag_coeff) + "\n"
    
    # FieldData
    fd_el = ET.SubElement(ug, "FieldData")
    st_arr = ET.SubElement(fd_el, "DataArray", type="String", Name="solver_type", format="ascii")
    st_arr.text = f"\n{solver_type}\n"
    
    if coefficients is not None:
        for name, val in coefficients.items():
            arr = ET.SubElement(fd_el, "DataArray", type="Float32", Name=name, format="ascii")
            arr.text = f"\n{val:.6f}\n"
            
    tree = ET.ElementTree(vtk_file)
    tree.write(filepath, encoding="utf-8", xml_declaration=True)

def read_vtu(filepath):
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    ug = root.find("UnstructuredGrid")
    piece = ug.find("Piece")
    
    n_points = int(piece.attrib["NumberOfPoints"])
    n_cells = int(piece.attrib["NumberOfCells"])
    
    # Points
    pts_el = piece.find("Points").find("DataArray")
    points = np.fromstring(pts_el.text, sep=" ").reshape((n_points, 3))
    
    # Cells
    cells_el = piece.find("Cells")
    connectivity = []
    offsets = []
    types = []
    for arr in cells_el.findall("DataArray"):
        name = arr.attrib.get("Name")
        if name == "connectivity":
            connectivity = np.fromstring(arr.text, dtype=int, sep=" ")
        elif name == "offsets":
            offsets = np.fromstring(arr.text, dtype=int, sep=" ")
        elif name == "types":
            types = np.fromstring(arr.text, dtype=int, sep=" ")
            
    # CellData
    cd_el = piece.find("CellData")
    cell_data = {}
    for arr in cd_el.findall("DataArray"):
        name = arr.attrib.get("Name")
        num_components = int(arr.attrib.get("NumberOfComponents", 1))
        data = np.fromstring(arr.text, sep=" ")
        if num_components > 1:
            data = data.reshape((n_cells, num_components))
        cell_data[name] = data
        
    # FieldData
    field_data = {}
    fd_el = ug.find("FieldData")
    if fd_el is not None:
        for arr in fd_el.findall("DataArray"):
            name = arr.attrib.get("Name")
            field_data[name] = arr.text.strip()
            
    return points, connectivity, offsets, types, cell_data, field_data

def run_post_process(analysis_folder):
    analysis_folder = Path(analysis_folder)
    vlm_vtu = analysis_folder / "results_vlm.vtu"
    nvlm_vtu = analysis_folder / "results_nvlm.vtu"
    
    vtu_files = []
    if nvlm_vtu.is_file():
        vtu_files.append((nvlm_vtu, "NL-VLM Solver Results"))
    if vlm_vtu.is_file():
        vtu_files.append((vlm_vtu, "VLM Solver Results"))
        
    if not vtu_files:
        print(f"{RED}Error: No results_vlm.vtu or results_nvlm.vtu found in {analysis_folder}. Run a solver first.{RESET}")
        return
        
    if len(vtu_files) == 1:
        vtu_path, desc = vtu_files[0]
    else:
        print(f"\n{BOLD}Select VTU file to post-process:{RESET}")
        for idx, (_, desc) in enumerate(vtu_files):
            print(f"  [{CYAN}{idx}{RESET}] {desc}")
        try:
            choice = int(input(f"Select choice (0-1): ").strip())
            if 0 <= choice < len(vtu_files):
                vtu_path, _ = vtu_files[choice]
            else:
                print(f"{RED}Invalid index. Using latest results.{RESET}")
                vtu_path, _ = vtu_files[0]
        except ValueError:
            vtu_path, _ = vtu_files[0]

    print(f"\n{YELLOW}Reading panel results from {vtu_path} and generating plots...{RESET}")
    points, connectivity, offsets, types, cell_data, field_data = read_vtu(vtu_path)
    solver_type = field_data.get("solver_type", "VLM").upper().replace("-", "")
    
    # Read dynamics parameters from FieldData
    v_inf = float(field_data.get("v_inf", 1.0))
    alfa_deg = float(field_data.get("alfa_deg", 5.0))
    
    setup_plot_style()
    plots_dir = analysis_folder / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    n_cells = len(types)
    gammas = cell_data["gamma"]
    
    # 1. 3D Surface Heatmap
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection='3d')
    title_suffix = "(NL-VLM)" if solver_type == "NLVLM" else "(VLM)"
    ax.set_title(f'3D Wing Surface Mesh colored by Vortex Strength {title_suffix}', 
                 fontsize=13, fontweight='bold', pad=20)
    
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    ax.xaxis._axinfo["grid"]['color'] = (0.8, 0.8, 0.8, 0.4)
    ax.yaxis._axinfo["grid"]['color'] = (0.8, 0.8, 0.8, 0.4)
    ax.zaxis._axinfo["grid"]['color'] = (0.8, 0.8, 0.8, 0.4)
    
    ax.set_xlabel('x [m]', fontsize=10, labelpad=10)
    ax.set_ylabel('y [m]', fontsize=10, labelpad=10)
    ax.set_zlabel('z [m]', fontsize=10, labelpad=10)
    
    norm = colors.Normalize(vmin=min(gammas), vmax=max(gammas))
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm)
    
    poly_list = []
    face_colors = []
    for k in range(n_cells):
        pts_indices = connectivity[4*k : 4*k+4]
        poly = points[pts_indices].tolist()
        poly_list.append(poly)
        face_colors.append(mapper.to_rgba(gammas[k]))
        
    mesh = Poly3DCollection(poly_list, edgecolor=(0.2, 0.2, 0.2, 0.2), lw=0.4, alpha=0.9)
    mesh.set_facecolor(face_colors)
    ax.add_collection3d(mesh)
    
    x_coords, y_coords, z_coords = points[:, 0], points[:, 1], points[:, 2]
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_mid = 0.5 * (max(z_coords) + min(z_coords))
    z_lim_half = max((max(z_coords) - min(z_coords)) * 0.5, 0.05 * (y_max - y_min))
    ax.set_zlim(z_mid - z_lim_half, z_mid + z_lim_half)
    
    ax.view_init(elev=20, azim=-125)
    ax.set_box_aspect((x_max - x_min, y_max - y_min, max(max(z_coords) - min(z_coords), 0.15 * (y_max - y_min))))
    
    mapper.set_array(gammas)
    cbar = fig.colorbar(mapper, ax=ax, shrink=0.55, aspect=18, pad=0.08)
    cbar.set_label(r'Vortex Strength $\Gamma$ [$m^2/s$]', fontsize=10, labelpad=10)
    
    heatmap_prefix = 'nlvlm_3d_heatmap' if solver_type == "NLVLM" else 'vlm_3d_heatmap'
    fig.savefig(plots_dir / f"{heatmap_prefix}.pdf", bbox_inches='tight')
    fig.savefig(plots_dir / f"{heatmap_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # 2. Spanwise Circulation
    y_centers = [np.mean(points[connectivity[4*k : 4*k+4], 1]) for k in range(n_cells)]
    gamma_dict = {}
    for y_c, g in zip(y_centers, gammas):
        key = round(y_c, 5)
        gamma_dict[key] = gamma_dict.get(key, []) + [g]
            
    sorted_keys = sorted(gamma_dict.keys())
    le_gamma = [sum(gamma_dict[k]) for k in sorted_keys]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if solver_type == "NLVLM":
        ax.plot(sorted_keys, le_gamma, 'g-o', lw=1.8, ms=4, label='NL-VLM')
        ax.set_title('Spanwise Circulation Distribution (NL-VLM)')
        circ_prefix = 'nlvlm_circulation'
    else:
        ax.plot(sorted_keys, le_gamma, 'r-o', lw=1.8, ms=4, label='VLM')
        ax.set_title('Leading Edge Total Circulation Distribution (VLM)')
        circ_prefix = 'le_circulation_3d_vlm'
        
    ax.set_xlabel('Spanwise position y [m]')
    ax.set_ylabel('gamma [m^2/s]')
    ax.grid(True, linestyle='--', alpha=0.3)
    fig.savefig(plots_dir / f"{circ_prefix}.pdf", bbox_inches='tight')
    fig.savefig(plots_dir / f"{circ_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Reconstruct section grid to calculate spanwise Cl, Cdi, Cdp
    first_y = y_centers[0]
    n_chord_sections = sum(1 for y in y_centers if abs(y - first_y) < 1e-3)
    n_span_sections = n_cells // n_chord_sections
    
    y_dist = []
    cl_dist = []
    cdi_dist = []
    cdp_dist = []
    
    alpha_rad = np.radians(alfa_deg)
    cos_a, sin_a = np.cos(alpha_rad), np.sin(alpha_rad)
    
    for i in range(n_span_sections):
        indices = range(i * n_chord_sections, (i + 1) * n_chord_sections)
        y_c = np.mean([y_centers[idx] for idx in indices])
        y_dist.append(y_c)
        
        S_i = sum(cell_data["area"][idx] for idx in indices)
        F_body = np.sum([cell_data["force"][idx] for idx in indices], axis=0)
        
        # Rotate body force to stability axes
        L_i = -F_body[0] * sin_a + F_body[2] * cos_a
        D_induced_i = F_body[0] * cos_a + F_body[2] * sin_a
        
        # Normalize by 0.5 * rho * Vinf^2 * S_i (with rho=1.0)
        cl_i = (2.0 * L_i) / (v_inf**2 * S_i)
        cdi_i = (2.0 * D_induced_i) / (v_inf**2 * S_i)
        cdp_i = np.mean([cell_data.get("profile_drag_coeff", np.zeros(n_cells))[idx] for idx in indices])
        
        cl_dist.append(cl_i)
        cdi_dist.append(cdi_i)
        cdp_dist.append(cdp_i)
        
    y_dist = np.array(y_dist)
    cl_dist = np.array(cl_dist)
    cdi_dist = np.array(cdi_dist)
    cdp_dist = np.array(cdp_dist)
    
    # 3. Plot Spanwise Local Lift Coefficient Cl
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(y_dist, cl_dist, 'b-d', lw=1.8, ms=4, label=f'Local $C_l$ ({solver_type})')
    ax.set_title(f'Spanwise Local Lift Coefficient Distribution {title_suffix}')
    ax.set_xlabel('Spanwise position y [m]')
    ax.set_ylabel('Sectional Lift Coefficient $C_l$')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend()
    prefix = 'nvlm' if solver_type == "NLVLM" else 'vlm'
    fig.savefig(plots_dir / f"{prefix}_spanwise_lift_distribution.pdf", bbox_inches='tight')
    fig.savefig(plots_dir / f"{prefix}_spanwise_lift_distribution.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # 4. Plot Spanwise Local Drag Coefficients Cdi and Cdp
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(y_dist, cdi_dist, 'r-s', lw=1.5, ms=4, label='Induced $C_{di}$')
    if solver_type == "NLVLM":
        ax.plot(y_dist, cdp_dist, 'g-^', lw=1.5, ms=4, label='Profile $C_{dp}$')
        ax.plot(y_dist, cdi_dist + cdp_dist, 'k--', lw=1.8, label='Total Section $C_d$')
    ax.set_title(f'Spanwise Local Drag Coefficient Distribution {title_suffix}')
    ax.set_xlabel('Spanwise position y [m]')
    ax.set_ylabel('Sectional Drag Coefficient $C_d$')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend()
    fig.savefig(plots_dir / f"{prefix}_spanwise_drag_distribution.pdf", bbox_inches='tight')
    fig.savefig(plots_dir / f"{prefix}_spanwise_drag_distribution.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"{GREEN}Post-processing complete. Output plots saved to plots/ directory.{RESET}")

def run_plot_param(config, analysis_folder):
    print(f"\n{YELLOW}Generating wing parameterization diagram...{RESET}")
    stations = sorted(config["wing"]["stations"], key=lambda st: float(st["y"]))
    
    setup_plot_style()
    plots_dir = Path(analysis_folder) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    fig = plt.figure(figsize=(12, 10))
    
    # 1. PLANFORM VIEW
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax1.set_title("Wing Planform View (Top-Down Parameterization)", fontsize=13, fontweight='bold', pad=15)
    
    y_stations = np.array([st["y"] for st in stations])
    x_le = np.array([st["x_le"] for st in stations])
    chords = np.array([st["chord"] for st in stations])
    x_te = x_le + chords
    
    y_full = np.concatenate([-y_stations[::-1], y_stations])
    x_le_full = np.concatenate([x_le[::-1], x_le])
    x_te_full = np.concatenate([x_te[::-1], x_te])
    
    ax1.plot(y_full, x_le_full, 'k-', lw=2, label="Leading Edge (LE)")
    ax1.plot(y_full, x_te_full, 'k--', lw=1.5, label="Trailing Edge (TE)")
    ax1.plot([y_full[0], y_full[0]], [x_le_full[0], x_te_full[0]], 'k-', lw=2)
    ax1.plot([y_full[-1], y_full[-1]], [x_le_full[-1], x_te_full[-1]], 'k-', lw=2)
    ax1.axvline(0, color='gray', linestyle=':', alpha=0.7, label="Symmetry Line")
    
    for idx, st in enumerate(stations):
        y_val, x_start, c_val = st["y"], st["x_le"], st["chord"]
        ax1.plot([y_val, y_val], [x_start, x_start + c_val], 'r-', alpha=0.8, lw=1.5)
        ax1.annotate(f"$c_{idx}$ = {c_val}m\n(Station {idx})", 
                     xy=(y_val, x_start + 0.5 * c_val), 
                     xytext=(y_val + 0.25, x_start + 0.5 * c_val),
                     arrowprops=dict(arrowstyle="->", color="red", lw=0.8),
                     fontsize=9, color="darkred")
        ax1.axvline(y_val, color='blue', linestyle='-.', alpha=0.3, lw=1)
        if idx > 0:
            ax1.text(y_val, x_le_full.min() - 0.2, f"$y_{idx}$ = {y_val}m", ha='center', va='top', fontsize=9, color='blue')
            ax1.axvline(-y_val, color='blue', linestyle='-.', alpha=0.15, lw=1)
            
    ax1.text(0.0, x_le_full.min() - 0.2, "$y_0$ = 0.0m\n(Root)", ha='center', va='top', fontsize=9, color='blue')
    span = 2.0 * y_stations[-1]
    ax1.annotate('', xy=(-span/2, x_te_full.max() + 0.3), xytext=(span/2, x_te_full.max() + 0.3),
                 arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax1.text(0.0, x_te_full.max() + 0.4, f"Total Wingspan $b$ = {span}m", ha='center', va='bottom', fontweight='bold')
    
    ax1.set_xlabel("Spanwise location y [m]")
    ax1.set_ylabel("Chordwise location x [m]")
    ax1.invert_yaxis()
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper left")
    
    # 2. DIHEDRAL VIEW
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    ax2.set_title("Front View (Dihedral & Station Elevations)", fontsize=11, fontweight='bold', pad=10)
    z_le = np.array([st.get("z_le", 0.0) for st in stations])
    y_full_di = np.concatenate([-y_stations[::-1], y_stations])
    z_le_full = np.concatenate([z_le[::-1], z_le])
    
    ax2.plot(y_full_di, z_le_full, 'b-o', ms=4, label="LE Elevation ($z_{le}$)")
    for idx, st in enumerate(stations):
        y_val, z_val = st["y"], st.get("z_le", 0.0)
        if idx > 0:
            ax2.text(y_val, z_val + 0.05, f"$z_{idx}$={z_val}m", ha='center', va='bottom', fontsize=8, color='blue')
            ax2.text(-y_val, z_val + 0.05, f"$z_{idx}$={z_val}m", ha='center', va='bottom', fontsize=8, color='blue', alpha=0.6)
    ax2.text(0.0, z_le[0] + 0.05, f"$z_0$={z_le[0]}m", ha='center', va='bottom', fontsize=8, color='blue')
    ax2.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Spanwise location y [m]")
    ax2.set_ylabel("Elevation z [m]")
    ax2.grid(True, alpha=0.3)
    
    # 3. TWIST SCHEMA
    ax3 = plt.subplot2grid((2, 2), (1, 1))
    ax3.set_title("Station Section Parameters (Airfoil & Twist)", fontsize=11, fontweight='bold', pad=10)
    
    theta_deg = 15.0
    theta_rad = np.radians(theta_deg)
    c_illustration = 1.0
    x_chord_line = np.array([0.0, c_illustration * np.cos(-theta_rad)])
    z_chord_line = np.array([0.0, c_illustration * np.sin(-theta_rad)])
    
    ax3.plot([0.0, c_illustration + 0.2], [0.0, 0.0], 'k:', alpha=0.6, label="Horizontal Ref")
    ax3.plot(x_chord_line, z_chord_line, 'r-', lw=2, label="Twisted Chord Line")
    
    angle_arc_x = 0.4 * np.cos(np.linspace(0, -theta_rad, 20))
    angle_arc_z = 0.4 * np.sin(np.linspace(0, -theta_rad, 20))
    ax3.plot(angle_arc_x, angle_arc_z, 'k-', lw=1.0)
    ax3.text(0.45, -0.08, r"Twist $\theta_{deg}$", ha='left', va='top', fontweight='bold')
    
    s_camber = np.linspace(0.0, 1.0, 50)
    z_camber_val = 0.08 * np.sin(np.pi * s_camber)
    x_camber_rot = s_camber * c_illustration * np.cos(-theta_rad) - z_camber_val * np.sin(-theta_rad)
    z_camber_rot = s_camber * c_illustration * np.sin(-theta_rad) + z_camber_val * np.cos(-theta_rad)
    ax3.plot(x_camber_rot, z_camber_rot, 'g--', lw=1.5, label="Camber Line")
    
    ax3.text(0.0, 0.05, "LE ($x_{le}, y, z_{le}$)", ha='right', va='bottom', fontsize=9, fontweight='bold')
    ax3.text(x_chord_line[1], z_chord_line[1] - 0.05, "TE", ha='left', va='top', fontsize=9, fontweight='bold')
    ax3.text(0.5 * x_chord_line[1], 0.5 * z_chord_line[1] - 0.15, "Chord length $c$", ha='center', va='top', color='red')
    
    text_y_offset = -0.5
    ax3.text(0.0, text_y_offset, "Actual Station Twist & Airfoils:", fontsize=9, fontweight='bold')
    for idx, st in enumerate(stations):
        tw, af = st.get("twist_deg", 0.0), st.get("airfoil_name", "None")
        ax3.text(0.0, text_y_offset - 0.08 * (idx + 1), f"  Station {idx}: $\\theta$ = {tw}°, Airfoil: {af}", fontsize=8.5)
        
    ax3.set_xlim(-0.2, c_illustration + 0.3)
    ax3.set_ylim(text_y_offset - 0.08 * (len(stations) + 2), 0.3)
    ax3.axis('off')
    ax3.legend(loc="upper right")
    
    fig.tight_layout()
    pdf_path = plots_dir / "wing_parameterization.pdf"
    png_path = plots_dir / "wing_parameterization.png"
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"{GREEN}Success! Wing parameterization diagram saved to {BOLD}{pdf_path}{RESET}")

def print_ascii_planform(mesher):
    width = 50
    height = 15
    grid_chars = [[" " for _ in range(width)] for _ in range(height)]
    
    span_half = mesher.span / 2.0
    
    x_coords = []
    for col in range(width):
        y = col * (span_half / (width - 1))
        chord, x_le, z_le, twist_rad, z_camber = mesher._interpolate_station_properties(y, 0.0)
        x_coords.append(x_le)
        x_coords.append(x_le + chord)
    x_min, x_max = min(x_coords), max(x_coords)
    x_range = x_max - x_min if x_max - x_min > 1e-6 else 1.0
    
    y_nodes = np.linspace(0.0, span_half, mesher.n_span)
    s_nodes = np.linspace(0.0, 1.0, mesher.n_chord)
    
    # Draw horizontal chordwise lines (ribs)
    for y_n in y_nodes:
        col = int(round(y_n * (width - 1) / span_half))
        col = max(0, min(width - 1, col))
        
        chord, x_le, z_le, twist_rad, z_camber = mesher._interpolate_station_properties(y_n, 0.0)
        for s in np.linspace(0.0, 1.0, 50):
            x = x_le + s * chord
            row = int(round((x - x_min) * (height - 1) / x_range))
            row = max(0, min(height - 1, row))
            grid_chars[row][col] = "|"
            
    # Draw spanwise panel lines
    for s_n in s_nodes:
        for col in range(width):
            y = col * (span_half / (width - 1))
            chord, x_le, z_le, twist_rad, z_camber = mesher._interpolate_station_properties(y, s_n)
            row = int(round((x_le - x_min) * (height - 1) / x_range))
            row = max(0, min(height - 1, row))
            
            if s_n == 0.0:
                grid_chars[row][col] = "=" 
            elif s_n == 1.0:
                grid_chars[row][col] = "=" 
            else:
                if grid_chars[row][col] == " ":
                    grid_chars[row][col] = "-"
                    
    print(f"\n{BOLD}{CYAN}--- 2D HALF-WING PLANFORM GRID DISCRETIZATION ---{RESET}")
    print("Chordwise (x) pointing Down, Spanwise (y) pointing Right")
    print("LE (Top boundary, '='), TE (Bottom boundary, '='), Ribs ('|'), Panel Lines ('-')\n")
    print("Root (y=0)" + " " * (width - 18) + "Tip (y=span/2)")
    print("   +" + "-" * width + "+")
    for r in range(height):
        row_str = "".join(grid_chars[r])
        print(f" {r:2d} |{row_str}|")
    print("   +" + "-" * width + "+")

def print_ascii_3d_mean_line(mesher):
    width = 55
    height = 15
    grid_chars = [[" " for _ in range(width)] for _ in range(height)]
    
    y_vals = np.linspace(-mesher.span / 2.0, mesher.span / 2.0, 60)
    
    points_3d = []
    for y in y_vals:
        chord, x_le, z_le, twist_rad, z_camber = mesher._interpolate_station_properties(y, 0.25)
        dx = 0.25 * chord
        dz = z_camber
        x_rot = dx * np.cos(twist_rad) - dz * np.sin(twist_rad)
        z_rot = dx * np.sin(twist_rad) + dz * np.cos(twist_rad)
        points_3d.append([x_le + x_rot, y, z_le + z_rot])
        
    points_3d = np.array(points_3d)
    
    cos30 = np.cos(np.radians(30))
    sin30 = np.sin(np.radians(30))
    
    u_proj = points_3d[:, 1] * cos30 - points_3d[:, 0] * sin30
    v_proj = points_3d[:, 2] + points_3d[:, 1] * sin30 + points_3d[:, 0] * cos30
    
    u_min, u_max = min(u_proj), max(u_proj)
    v_min, v_max = min(v_proj), max(v_proj)
    
    u_range = u_max - u_min if u_max - u_min > 1e-6 else 1.0
    v_range = v_max - v_min if v_max - v_min > 1e-6 else 1.0
    
    for u, v in zip(u_proj, v_proj):
        col = int(round((u - u_min) * (width - 1) / u_range))
        row = int(round((v - v_min) * (height - 1) / v_range))
        
        col = max(0, min(width - 1, col))
        row = max(0, min(height - 1, row))
        
        grid_chars[(height - 1) - row][col] = "o"
        
    left_tip_col = int(round((u_proj[0] - u_min) * (width - 1) / u_range))
    left_tip_row = (height - 1) - int(round((v_proj[0] - v_min) * (height - 1) / v_range))
    right_tip_col = int(round((u_proj[-1] - u_min) * (width - 1) / u_range))
    right_tip_row = (height - 1) - int(round((v_proj[-1] - v_min) * (height - 1) / v_range))
    root_col = int(round((u_proj[30] - u_min) * (width - 1) / u_range))
    root_row = (height - 1) - int(round((v_proj[30] - v_min) * (height - 1) / v_range))
    
    grid_chars[left_tip_row][left_tip_col] = "L"
    grid_chars[right_tip_row][right_tip_col] = "R"
    grid_chars[root_row][root_col] = "C"
    
    print(f"\n{BOLD}{CYAN}--- 3D MEAN CAMBER LINE ISOMETRIC PROJECTION ---{RESET}")
    print("Center/Root ('C'), Left Tip ('L'), Right Tip ('R'), Path ('o')\n")
    print("   +" + "-" * width + "+")
    for r in range(height):
        row_str = "".join(grid_chars[r])
        print(f"   |{row_str}|")
    print("   +" + "-" * width + "+")

def show_ascii_plots(config):
    mesher = WingMesher(config, data_folder="airfoilDATA")
    print_ascii_planform(mesher)
    print_ascii_3d_mean_line(mesher)

def modify_config_live(config):
    while True:
        print(f"\n{BOLD}Live Configuration Submenu (In-Memory Only):{RESET}")
        print(f"  [{CYAN}1{RESET}] Change Angle of Attack (alpha) [current: {config['solver']['alfa_deg']}°]")
        print(f"  [{CYAN}2{RESET}] Change Freestream Velocity (Vinf) [current: {config['solver']['v_inf']} m/s]")
        print(f"  [{CYAN}3{RESET}] Change Mesh Grid Density [current: {config['meshing']['n_span']} span x {config['meshing']['n_chord']} chord]")
        print(f"  [{CYAN}4{RESET}] Change Station Twist/Chord")
        print(f"  [{CYAN}5{RESET}] Change NL-VLM Convergence Tolerance [current: {config['solver'].get('tol', 1e-3)}]")
        print(f"  [{CYAN}6{RESET}] Back to Main Menu")
        
        sub_choice = input(f"\nSelect an option (1-6): ").strip()
        if sub_choice == "1":
            try:
                new_alpha = float(input(f"Enter new Angle of Attack in degrees: ").strip())
                config["solver"]["alfa_deg"] = new_alpha
                print(f"{GREEN}Angle of Attack updated to {new_alpha}°{RESET}")
            except ValueError:
                print(f"{RED}Invalid input. Please enter a number.{RESET}")
        elif sub_choice == "2":
            try:
                new_vinf = float(input(f"Enter new Freestream Velocity: ").strip())
                config["solver"]["v_inf"] = new_vinf
                print(f"{GREEN}Freestream velocity updated to {new_vinf} m/s{RESET}")
            except ValueError:
                print(f"{RED}Invalid input. Please enter a number.{RESET}")
        elif sub_choice == "3":
            try:
                new_span = int(input(f"Enter new number of spanwise nodes (n_span): ").strip())
                new_chord = int(input(f"Enter new number of chordwise nodes (n_chord): ").strip())
                config["meshing"]["n_span"] = new_span
                config["meshing"]["n_chord"] = new_chord
                print(f"{GREEN}Mesh density updated to {new_span} x {new_chord}{RESET}")
            except ValueError:
                print(f"{RED}Invalid input. Please enter integers.{RESET}")
        elif sub_choice == "4":
            print(f"\n{BOLD}Select wing station to modify:{RESET}")
            stations = config["wing"]["stations"]
            for idx, st in enumerate(stations):
                print(f"  Station {idx}: y={st['y']}m, chord={st['chord']}m, twist={st['twist_deg']}°, airfoil={st.get('airfoil_name', 'None')}")
            try:
                st_idx = int(input(f"Enter station index: ").strip())
                if 0 <= st_idx < len(stations):
                    st = stations[st_idx]
                    chord_input = input(f"Enter new chord [{st['chord']}m] (leave empty to keep): ").strip()
                    twist_input = input(f"Enter new twist [{st['twist_deg']}°] (leave empty to keep): ").strip()
                    
                    if chord_input:
                        st["chord"] = float(chord_input)
                    if twist_input:
                        st["twist_deg"] = float(twist_input)
                    print(f"{GREEN}Station {st_idx} updated!{RESET}")
                else:
                    print(f"{RED}Invalid index.{RESET}")
            except ValueError:
                print(f"{RED}Invalid input.{RESET}")
        elif sub_choice == "5":
            try:
                new_tol = float(input(f"Enter new NL-VLM Convergence Tolerance (e.g. 1e-4): ").strip())
                config["solver"]["tol"] = new_tol
                print(f"{GREEN}Tolerance updated to {new_tol}{RESET}")
            except ValueError:
                print(f"{RED}Invalid input. Please enter a float.{RESET}")
        elif sub_choice == "6":
            break

def run_section_analysis(config, analysis_folder):
    print(f"\n{BOLD}Wing Stations for Polar Analysis:{RESET}")
    stations = config["wing"]["stations"]
    for idx, st in enumerate(stations):
        print(f"  [{CYAN}{idx}{RESET}] Station {idx}: y={st['y']}m, chord={st['chord']}m, airfoil={st.get('airfoil_name', 'None')}")
    print(f"  [{CYAN}c{RESET}] Custom Section Analysis (specify custom parameters)")
    print(f"  [{CYAN}e{RESET}] Exit / Back to Main Menu")
        
    choice = input(f"\nSelect a station index, 'c' or 'e': ").strip().lower()
    
    if choice in ("e", "exit"):
        return
    elif choice in ("c", "custom"):
        af_name = input("Enter custom airfoil name (e.g. naca2412, ag18): ").strip().lower()
        if not af_name:
            print(f"{RED}Airfoil name cannot be empty.{RESET}")
            return
        
        chord_str = input("Enter chord in meters [default: 1.0 m]: ").strip()
        chord_val = float(chord_str) if chord_str else 1.0
        
        vinf_str = input(f"Enter freestream velocity in m/s [default: {config['solver']['v_inf']} m/s]: ").strip()
        v_inf = float(vinf_str) if vinf_str else float(config["solver"]["v_inf"])
        st_idx = "custom"
    else:
        try:
            st_idx = int(choice)
            if not (0 <= st_idx < len(stations)):
                print(f"{RED}Invalid index.{RESET}")
                return
        except ValueError:
            print(f"{RED}Invalid input.{RESET}")
            return
            
        st = stations[st_idx]
        af_name = st.get("airfoil_name")
        if not af_name or af_name.lower() == "none":
            print(f"{RED}This station does not have an airfoil associated with it.{RESET}")
            return
            
        chord_val = float(st["chord"])
        
        # Ask for Vinf (default to config value)
        vinf_str = input(f"Enter freestream velocity in m/s [default: {config['solver']['v_inf']} m/s]: ").strip()
        v_inf = float(vinf_str) if vinf_str else float(config["solver"]["v_inf"])
    
    re = (v_inf * chord_val) / 1.46e-5
    re_target = min([50000, 100000, 200000, 500000, 1000000], key=lambda x: abs(x - re))
    
    print(f"\n{GREEN}Analyzing Station {st_idx} ({af_name}) at Vinf = {v_inf:.2f} m/s:{RESET}")
    print(f"  Local Chord = {chord_val:.2f}m")
    print(f"  Calculated Re = {re:,.0f} -> Matched Database Re = {re_target:,}")
    
    # Load AeroFoil polar
    af = AeroFoil(af_name, "airfoilDATA")
    af.downloadPerformance(re_target)
    
    cl_fit_func = af.fitCL(5)
    cd_fit_func = af.fitCD(5)
    
    plots_dir = Path(analysis_folder) / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    while True:
        print(f"\n{BOLD}Section Analysis Options:{RESET}")
        print(f"  [{CYAN}1{RESET}] Export Fitted Polar Data to CSV")
        print(f"  [{CYAN}2{RESET}] Plot Fitted Polar Curves (PDF & PNG)")
        print(f"  [{CYAN}3{RESET}] Plot Airfoil Camber Line (PDF & PNG)")
        print(f"  [{CYAN}4{RESET}] Return to Main Menu")
        
        sub_choice = input(f"\nSelect an option (1-4): ").strip()
        
        if sub_choice == "1":
            # Generate polar data points from -10 to 15 degrees
            aoas = np.linspace(-10.0, 15.0, 51)
            csv_path = plots_dir / f"station_{st_idx}_polar_{re_target}.csv"
            
            with open(csv_path, "w") as f:
                f.write("AoA_deg,Cl_fitted,Cd_fitted\n")
                for aoa in aoas:
                    cl_val = cl_fit_func(radians(aoa))
                    cd_val = cd_fit_func(radians(aoa))
                    f.write(f"{aoa:.2f},{cl_val:.6f},{cd_val:.6f}\n")
            print(f"{GREEN}Fitted polar data successfully exported to {BOLD}{csv_path}{RESET}")
            
        elif sub_choice == "2":
            setup_plot_style()
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
            
            # AoA range for plotting fit
            aoas_fit_deg = np.linspace(min(af.alfa), max(af.alfa), 120)
            cl_fit_vals = [cl_fit_func(radians(a)) for a in aoas_fit_deg]
            cd_fit_vals = [cd_fit_func(radians(a)) for a in aoas_fit_deg]
            
            # Left panel: Cl vs alpha
            ax1.scatter(af.alfa, af.CL, color='gray', alpha=0.6, s=15, label='Raw Database Polar')
            ax1.plot(aoas_fit_deg, cl_fit_vals, 'b-', lw=2, label='Cubic Spline Fit')
            ax1.set_title(f'Lift Coefficient $C_l$ vs $\\alpha$ ({af_name.upper()} @ Re={re_target:,})', fontsize=11, fontweight='bold')
            ax1.set_xlabel('Angle of Attack $\\alpha$ [deg]')
            ax1.set_ylabel('Lift Coefficient $C_l$')
            ax1.grid(True, linestyle='--', alpha=0.3)
            ax1.legend()
            
            # Right panel: Cd vs alpha
            ax2.scatter(af.alfa, af.CD, color='gray', alpha=0.6, s=15, label='Raw Database Polar')
            ax2.plot(aoas_fit_deg, cd_fit_vals, 'r-', lw=2, label='Cubic Spline Fit')
            ax2.set_title(f'Drag Coefficient $C_d$ vs $\\alpha$ ({af_name.upper()} @ Re={re_target:,})', fontsize=11, fontweight='bold')
            ax2.set_xlabel('Angle of Attack $\\alpha$ [deg]')
            ax2.set_ylabel('Drag Coefficient $C_d$')
            ax2.grid(True, linestyle='--', alpha=0.3)
            ax2.legend()
            
            pdf_path = plots_dir / f"station_{st_idx}_polar_{re_target}.pdf"
            png_path = plots_dir / f"station_{st_idx}_polar_{re_target}.png"
            fig.savefig(pdf_path, bbox_inches='tight')
            fig.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            print(f"{GREEN}Fitted polar curves plotted successfully!{RESET}")
            print(f"  PDF: {BOLD}{pdf_path}{RESET}")
            print(f"  PNG: {BOLD}{png_path}{RESET}")
            
        elif sub_choice == "3":
            setup_plot_style()
            fig, ax = plt.subplots(figsize=(9, 4.5))
            
            x_raw_unique = np.unique(af.x_base)
            camber_func = af.fitCamberLine()
            y_raw_unique = [camber_func(x) for x in x_raw_unique]
            
            # Spline fitted values for the continuous line
            s_fit = np.linspace(0.0, 1.0, 100)
            camber_fit = [camber_func(s) for s in s_fit]
            # -----------------------
            
            # Plot elements
            ax.plot(af.x_base, af.y_base, 'k-', lw=1.5, label=f'Airfoil Contour ({af_name.upper()})')
            ax.scatter(x_raw_unique, y_raw_unique, color='red', marker='x', s=20, zorder=3, label='Discrete Camber Midpoints')
            ax.plot(s_fit, camber_fit, 'g--', lw=2.2, label='Cubic Spline Camber Fit')
            
            ax.set_title(f'Mean Camber Line Fit and Airfoil Profile: {af_name.upper()}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Chordwise Fraction x/c')
            ax.set_ylabel('Thickness & Camber z/c')
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.axis('equal')
            ax.legend()

            pdf_path = plots_dir / f"station_{st_idx}_camber.pdf"
            png_path = plots_dir / f"station_{st_idx}_camber.png"
            fig.savefig(pdf_path, bbox_inches='tight')
            fig.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
        elif sub_choice == "4":
            break

def run_solver(solver_type, config, analysis_folder):
    print(f"\n{YELLOW}Starting {solver_type} analysis...{RESET}")
    mesher = WingMesher(config, data_folder="airfoilDATA")
    alfa_geom = radians(float(config["solver"]["alfa_deg"]))
    Vinf = float(config["solver"]["v_inf"])
    
    analysis_folder = Path(analysis_folder)
    analysis_folder.mkdir(parents=True, exist_ok=True)
    
    if solver_type == "VLM":
        vtu_path = analysis_folder / "results_vlm.vtu"
        print(f"{BLUE}[VLM]{RESET} Generating 3D panel mesh...")
        panels = mesher.generate_mesh()
        print(f"{BLUE}[VLM]{RESET} Mesh generated with {GREEN}{len(panels)}{RESET} panels.")

        print(f"{BLUE}[VLM]{RESET} Solving linear Vortex Lattice system...")
        vlm_solver = VortexLatticeSolver(panels)
        vlm_solver.setDynamicData(degrees(alfa_geom), Vinf)
        vlm_solver.solveSystem()

        vlm_CDi, vlm_CS, vlm_CL = vlm_solver.getAerodynamicCoefficients()
        
        # Normalize coefficients by Vinf^2 since solver doesn't normalize internally
        vlm_CDi = vlm_CDi / (Vinf**2)
        vlm_CL = vlm_CL / (Vinf**2)
        
        vlm_AR = (mesher.span**2) / vlm_solver.surfaceArea
        vlm_eff = (vlm_CL**2) / (np.pi * vlm_AR * vlm_CDi) if vlm_CDi > 1e-12 else 0.0

        print(f"""
{GREEN}----------------------------------------
VLM AERODYNAMIC COEFFICIENTS
  Surface Area  : {round(vlm_solver.surfaceArea, 4)} m^2
  Aspect Ratio  : {round(vlm_AR, 4)}
  Lift Coeff CL : {BOLD}{round(vlm_CL, 4)}{RESET}{GREEN}
  Drag Coeff CDi: {BOLD}{round(vlm_CDi, 4)}{RESET}{GREEN}
  Efficiency e  : {round(vlm_eff, 4)}
----------------------------------------{RESET}""")

        save_vtu(vtu_path, mesher, panels, "VLM", coefficients={
            "CL": vlm_CL,
            "CDi": vlm_CDi,
            "CDp": 0.0,
            "CD": vlm_CDi,
            "efficiency": vlm_eff,
            "v_inf": Vinf,
            "alfa_deg": float(config["solver"]["alfa_deg"])
        })
        print(f"{GREEN}Success! Panel data exported to {BOLD}{vtu_path}{RESET}")

    elif solver_type == "NLVLM":
        vtu_path = analysis_folder / "results_nvlm.vtu"
        print(f"{BLUE}[NL-VLM]{RESET} Generating 3D panel mesh...")
        panels = mesher.generate_mesh()
        
        damping = float(config["solver"].get("damping_factor", 0.05))
        tol = float(config["solver"].get("tol", 1e-3))
        max_iter = int(config["solver"].get("max_iter", 150))

        # Airfoil performance curves are dynamically downloaded/loaded based on local Reynolds numbers inside solve()

        print(f"{BLUE}[NL-VLM]{RESET} Solving coupled non-linear decambering loop...")
        nlvlm_solver = NonLinearVortexLatticeSolver(panels, mesher)
        control_pts, gamma_sec, iterations, error, nlvlm_CL, nlvlm_CDi, nlvlm_CDp, nlvlm_CD, nlvlm_eff = nlvlm_solver.solve(
            degrees(alfa_geom), Vinf, damping=damping, tol=tol, max_iter=max_iter
        )

        print(f"""
{GREEN}----------------------------------------
NL-VLM CONVERGENCE RESULTS
  Iterations        : {iterations}
  Final Residual    : {round(error, 5)}
----------------------------------------
NL-VLM AERODYNAMIC COEFFICIENTS
  Lift Coeff CL : {BOLD}{round(nlvlm_CL, 4)}{RESET}{GREEN}
  Induced CDi   : {BOLD}{round(nlvlm_CDi, 5)}{RESET}{GREEN}
  Profile CDp   : {BOLD}{round(nlvlm_CDp, 5)}{RESET}{GREEN}
  Total Drag CD : {BOLD}{round(nlvlm_CD, 5)}{RESET}{GREEN}
  Efficiency e  : {round(nlvlm_eff, 4)}
----------------------------------------{RESET}""")

        # Replicate station profile drag coefficients for individual panels in the VTU
        profile_drag_coeff = []
        for i in range(nlvlm_solver.n_span_sections):
            profile_drag_coeff.extend([nlvlm_solver.cd_profile_sec[i]] * nlvlm_solver.n_chord_sections)
            
        save_vtu(vtu_path, mesher, panels, "NLVLM", nlvlm_solver.decamber_angles, coefficients={
            "CL": nlvlm_CL,
            "CDi": nlvlm_CDi,
            "CDp": nlvlm_CDp,
            "CD": nlvlm_CD,
            "efficiency": nlvlm_eff,
            "v_inf": Vinf,
            "alfa_deg": float(config["solver"]["alfa_deg"])
        }, profile_drag_coeff=profile_drag_coeff)
        print(f"{GREEN}Success! Panel data exported to {BOLD}{vtu_path}{RESET}")

def load_config_from_path(config_path):
    p = Path(config_path)
    if not p.is_file():
        print(f"{RED}Error: Configuration file not found at {config_path}{RESET}")
        return None
    with open(p, "r") as f:
        return json.load(f)

def run_interactive_cli(config, analysis_folder):
    print(ASCII_LOGO)
    
    while True:
        print(f"\n{BOLD}Analysis Workspace Settings:{RESET}")
        print(f"  Analysis Folder: {CYAN}{analysis_folder}{RESET}")
        print(f"  Freestream AoA : {CYAN}{config['solver']['alfa_deg']}°{RESET} | Speed: {CYAN}{config['solver']['v_inf']} m/s{RESET}")
        
        print(f"\n{BOLD}Interactive Options:{RESET}")
        print(f"  [{CYAN}1{RESET}] Run Linear VLM Solver")
        print(f"  [{CYAN}2{RESET}] Run Non-Linear VLM Solver (NL-VLM)")
        print(f"  [{CYAN}3{RESET}] Run Post-Processing Plotter")
        print(f"  [{CYAN}4{RESET}] Generate Wing Parameterization PDF/PNG")
        print(f"  [{CYAN}5{RESET}] Show ASCII 2D Planform & 3D Mean Line Plots")
        print(f"  [{CYAN}6{RESET}] Modify Wing Configuration (Transient Session Editor)")
        print(f"  [{CYAN}7{RESET}] Perform Section Polar Analysis")
        print(f"  [{CYAN}8{RESET}] Exit")
        
        choice = input(f"\nSelect an option (1-8): ").strip()
        
        if choice == "1":
            run_solver("VLM", config, analysis_folder)
        elif choice == "2":
            run_solver("NLVLM", config, analysis_folder)
        elif choice == "3":
            run_post_process(analysis_folder)
        elif choice == "4":
            run_plot_param(config, analysis_folder)
        elif choice == "5":
            show_ascii_plots(config)
        elif choice == "6":
            modify_config_live(config)
        elif choice == "7":
            run_section_analysis(config, analysis_folder)
        elif choice == "8":
            print(f"\n{GREEN}Goodbye!{RESET}")
            break
        else:
            print(f"\n{RED}Invalid choice. Please input a number from 1 to 8.{RESET}")

def main():
    parser = argparse.ArgumentParser(description="Aerodynamic VLM & NL-VLM Solver Suite")
    parser.add_argument("-s", "--solver", choices=["VLM", "NL-VLM", "NLVLM"], help="Run specified solver directly")
    parser.add_argument("-p", "--postprocess", action="store_true", help="Run post-processing after solving")
    parser.add_argument("-c", "--config", help="Path to JSON configuration file")
    parser.add_argument("-a", "--analysis", help="Name of analysis folder to save results in")
    parser.add_argument("--param", action="store_true", help="Generate parameterization diagram")
    parser.add_argument("--ascii", action="store_true", help="Print ASCII 2D and 3D geometry plots")
    parser.add_argument("--cli", action="store_true", help="Force interactive CLI mode")
    
    args = parser.parse_args()
    
    # Prompt for configuration file path if not provided in arguments
    if args.config:
        config_path = args.config
    else:
        config_path_str = input(f"Enter configuration file path [default: config.json]: ").strip()
        config_path = config_path_str if config_path_str else "config.json"
        
    config = load_config_from_path(config_path)
    if config is None:
        sys.exit(1)
        
    # Prompt for analysis folder name if not provided in arguments
    if args.analysis:
        analysis_folder = args.analysis
    else:
        analysis_folder_str = input(f"Enter analysis folder name to save all results [default: analysis_results]: ").strip()
        analysis_folder = analysis_folder_str if analysis_folder_str else "analysis_results"
        
    analysis_folder = Path(analysis_folder)
    analysis_folder.mkdir(parents=True, exist_ok=True)
    
    if len(sys.argv) == 1 or args.cli or (not args.solver and not args.param and not args.ascii):
        run_interactive_cli(config, analysis_folder)
    else:
        print(ASCII_LOGO)
        if args.solver:
            solver_type = "NLVLM" if args.solver.upper() in ("NL-VLM", "NLVLM") else "VLM"
            run_solver(solver_type, config, analysis_folder)
            
        if args.postprocess:
            run_post_process(analysis_folder)
            
        if args.param:
            run_plot_param(config, analysis_folder)
            
        if args.ascii:
            show_ascii_plots(config)

if __name__ == "__main__":
    main()
