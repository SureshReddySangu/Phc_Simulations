import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import os 
import  sys

sys.path.append('/opt/lumerical/v231/api/python')
sys.path.append("opt/lumerical.v231.api")
sys.path.append("opt/lumerical/v231/bin")
import lumapi

#mirror parameters

H_steps:int= 5 
H_start:float= 500e-9 # in nm 
H_stop:float=  575e-9
H_values= np.linspace(H_start, H_stop,H_steps)
Amp_values:float = 220e-9
Z_thk:float  = 210e-9

open_cells:int = 5
unit_cells:int=20
close_cells:int = 7
total_cells:int = open_cells+unit_cells+close_cells
a:float = 320e-9
sample_dx:float = 1e-9
steps_per_cell :int = np.round(a/sample_dx)
n_Nitride:float= 1.99735

#strucutre  types
stru_type = ["cavity", "no_cavity"]
phc_type =  ["monogator", "holes", "fence_gate"]

#cavity parameters
cavity_rect_len:float = 90e-6 # in um 

#mirror_pads 
pad_rect_len:float = 2.5e-6

#mesh accuray 
mesh_acc:int = 2
fdtd_x_min= -pad_rect_len+0.5e-6 # in um 
margin:float = 2.5e-6

#file path 


file_path =  os.path.realpath(__file__)
dir_path = os.path.dirname(file_path)



# simlation check flag 
c :int = 299_792_458 #speed of the light
run_simulation:bool =1
simulation_time =  100*(2*cavity_rect_len)/c

wvl:float = 852e-9
bandwidth:float = 30e-9
low_wvl = wvl-bandwidth
high_wvl = wvl+bandwidth
num_points = (high_wvl-low_wvl)
freq_res= bandwidth*10 # using the 10 smapling poitns per wvl 
fdtd = lumapi.FDTD(hide=True)
#creat the structure

def make_phc_(phc, mod_rate, H):
    if phc == "monogator":
        print(mod_rate, "is mod_rate")
        x_pos = 0 # mirror starting point 
        total_vertices = int(total_cells*steps_per_cell)
        vertices = np.zeros((total_vertices,2))
        vertices_2 = np.zeros((total_vertices,2))
        vertex_index= 0
        for i in range(1, total_cells+1):
            if i <open_cells:
                modulation_rate = i*mod_rate/open_cells
            elif i<=open_cells+unit_cells:
                modulation_rate = mod_rate
            else:
                modulation_rate= mod_rate*(1-(i-(open_cells+unit_cells))/close_cells)
            for j in range(1, int(steps_per_cell+1)):
                x  = x_pos + j * sample_dx
                y_ = 0.5 * H + modulation_rate * np.sin(2 * np.pi * x / a) * 0.5
                vertices[vertex_index]   = (x,  y_)
                vertices_2[vertex_index] = (x, -y_)
                vertex_index += 1
            x_pos += a
        merged_polygon = np.vstack((vertices, np.flip(vertices_2, axis=0)))
        merged_polygon = np.flip(merged_polygon, axis=0)
        print(merged_polygon, "thisis")
        return merged_polygon
    elif phc=="fence_gate":
        print(mod_rate, "is mod_rate")
        x_pos = 0 # mirror starting point 
        total_vertices = int(total_cells*steps_per_cell)
        vertices = np.zeros((total_vertices,2))
        vertices_2 = np.zeros((total_vertices,2))
        vertex_index= 0
        for i in range(1, total_cells+1):
            if i <open_cells:
                modulation_rate = i*mod_rate/open_cells
            elif i<=open_cells+unit_cells:
                modulation_rate = mod_rate
            else:
                modulation_rate= mod_rate*(1-(i-(open_cells+unit_cells))/close_cells)
            for j in range(1, int(steps_per_cell+1)):
                x  = x_pos + j * sample_dx
                y_ = 0.5 * H + modulation_rate *np.sign(np.sin(2 * np.pi * x / a) * 0.5)
                vertices[vertex_index]   = (x,  y_)
                vertices_2[vertex_index] = (x, -y_)
                vertex_index += 1
            x_pos += a
        merged_polygon = np.vstack((vertices, np.flip(vertices_2, axis=0)))
        merged_polygon= np.flip(merged_polygon, axis=0)    
        print(merged_polygon, "thisis")
        return merged_polygon

# make the simualtion domain 
def make_fdtd_simu_domain(stru_type,phc, H,mod_rate,wg_thick, pad_rect_len, bare_mirror_len,):
    # we need to estimate the length of the device
    if stru_type =="no_cavity":
        if phc in ("monogator","fence_gate"):
            print("cavity_rect_len", cavity_rect_len)
            total_simu_domian_len = cavity_rect_len+bare_mirror_len*2+pad_rect_len*4
            print("total simultion domain size ", total_simu_domian_len)
            print(phc, "is Phc_type")
            merged_polygon = make_phc_(phc, mod_rate, H)
            # adding the rect
            fdtd.addrect() 
            fdtd.set("name", "left_mirror_left_pad") 
            fdtd.set("x min", -pad_rect_len)  # added -ve beaucse it is xmin 
            fdtd.set("x max", 0)   # since the mirror starts from the 0, so we extending the padding upto 0
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride) 
            #adding Monogator/ left mirror
            fdtd.addpoly() 
            fdtd.set("name", "Left_Mirror") 
            print(type(merged_polygon), "this is hte type of merged_polygon")
            fdtd.set("vertices", merged_polygon) 
            fdtd.set("z span", wg_thick) 
            fdtd.set("index", n_Nitride)

            #right padding rect
            fdtd.addrect() 
            fdtd.set("name", "left_mirror_right_pad") 
            fdtd.set("x min", bare_mirror_len) 
            fdtd.set("x max", bare_mirror_len+pad_rect_len)  # here we extended the mirror by pad_rect 
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride)


    else:
        if phc in ("monogator","fence_gate"):
            print("cavity_rect_len", cavity_rect_len)
            total_simu_domian_len = cavity_rect_len+bare_mirror_len*2+pad_rect_len*4
            print("total simultion domain size ", total_simu_domian_len)
            print(phc, "is Phc_type")
            merged_polygon = make_phc_(phc, mod_rate, H)
            # adding the rect
            fdtd.addrect() 
            fdtd.set("name", "left_mirror_left_pad") 
            fdtd.set("x min", -pad_rect_len)  # added -ve beaucse it is xmin 
            fdtd.set("x max", 0)   # since the mirror starts from the 0, so we extending the padding upto 0
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride) 
            #adding Monogator/ left mirror
            fdtd.addpoly() 
            fdtd.set("name", "Left_Mirror") 
            print(type(merged_polygon), "this is hte type of merged_polygon")
            fdtd.set("vertices", merged_polygon) 
            fdtd.set("z span", wg_thick) 
            fdtd.set("index", n_Nitride)

            #right padding rect
            fdtd.addrect() 
            fdtd.set("name", "left_mirror_right_pad") 
            fdtd.set("x min", bare_mirror_len) 
            fdtd.set("x max", bare_mirror_len+pad_rect_len)  # here we extended the mirror by pad_rect 
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride)

            #pad rect 
            fdtd.addrect()
            fdtd.set("name", "left_hole_mirror_right_padrect")
            fdtd.set("x min", total_cells*a+pad_rect_len)
            fdtd.set("x max", total_cells*a+2*pad_rect_len)
            fdtd.set("y min", -H/2)
            fdtd.set("y max", H/2)
            fdtd.set("z min", -wg_thick/2)
            fdtd.set("z max", wg_thick/2)
            fdtd.set("index", n_Nitride) 
            fdtd.select("left_hole_mirror_right_padrect")
            fdtd.addtogroup("miror_rect")
            #
            end_of_left_mirror = total_cells*a+2*pad_rect_len
            cavity_end = end_of_left_mirror+cavity_rect_len
            #
            #adding the cavity rec
            fdtd.addrect()
            fdtd.set("name", "cavity_rect")
            fdtd.set("x min", end_of_left_mirror)
            fdtd.set("x max",cavity_end)
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2)  
            fdtd.set("index", n_Nitride)

            fdtd.addrect() 
            fdtd.set("name", "right_mirror_left_pad") 
            fdtd.set("x min", cavity_end)  # added -ve beaucse it is xmin 
            fdtd.set("x max", cavity_end+pad_rect_len)   # since the mirror starts from the 0, so we extending the padding upto 0
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride) 
            #adding Monogator/ left mirror
            # Flip the polygon for the right mirror
            right_polygon = merged_polygon.copy()
            right_polygon[:, 0] = bare_mirror_len - right_polygon[:, 0]  # flip x about center

            fdtd.addpoly()
            fdtd.set("name", "Right_Mirror")
            fdtd.set("vertices", right_polygon)
            fdtd.set("z span", wg_thick)
            fdtd.set("index", n_Nitride)
            fdtd.set("x", cavity_end + pad_rect_len)
            end_of_structure = cavity_end+pad_rect_len+total_cells*a
           #right padding rect
            fdtd.addrect() 
            fdtd.set("name", "right_mirror_right_pad") 
            fdtd.set("x min",end_of_structure ) 
            fdtd.set("x max", end_of_structure+pad_rect_len)  # here we extended the mirror by pad_rect 
            fdtd.set("y min", -H/2) 
            fdtd.set("y max" , H/2) 
            fdtd.set( "z min", -wg_thick/2) 
            fdtd.set("z max", wg_thick/2) 
            fdtd.set("index", n_Nitride)





            #adding the fdtd
        source_ext_factor =6 #source extension factor 
        global source_x
        source_x = -pad_rect_len/2  #we place the source in the middle of the rectanlge 
        source_y_min = -source_ext_factor*(H)/2 
        source_y_max = source_ext_factor*(H)/2
        source_z_min = -source_ext_factor*wg_thick*0.5
        source_z_max = source_ext_factor*wg_thick*0.5 
        fdtd_y_min = -(H)-margin  
        fdtd_y_max = (H)+margin  

        fdtd_z_min = -(wg_thick) - margin 
        fdtd_z_max =  (wg_thick) + margin
        #R monitor
        global R_monitor_x
        R_monitor_x = source_x-0.5e-6  # we place the R monitor just before the fdtd. 
        R_monitor_y_min = fdtd_y_min 
        R_monitor_y_max = fdtd_y_max 
        R_monitor_z_min = fdtd_z_min 
        R_monitor_z_max = fdtd_z_max 
        # T monitor
        #we dynamically compute the T monitor x from the length of the strcure. 
        T_monitor_y_min = fdtd_y_min 
        T_monitor_y_max = fdtd_y_max 
        T_monitor_z_min = fdtd_z_min 
        T_monitor_z_max = fdtd_z_max 
        
        global fdtd_x_max
        if stru_type=="no_cavity":
            fdtd_x_max = total_simu_domian_len
        else:
            fdtd_x_max = end_of_structure+pad_rect_len
        #adding fdtd
        fdtd.addfdtd()
        fdtd.set("mesh accuracy", mesh_acc)
        fdtd.set("x min", fdtd_x_min) 
        fdtd.set("x max", fdtd_x_max)  # shoritng fdtd to cover the structuere fully 
        fdtd.set("y min", fdtd_y_min) 
        fdtd.set("y max", fdtd_y_max) 
        fdtd.set("z min", fdtd_z_min) 
        fdtd.set("z max", fdtd_z_max) 
        fdtd.set("use early shutoff", 0) 
        #adding the boundary codtions
        fdtd.set("y min bc", "Anti-Symmetric")
        fdtd.set("z min bc", "Symmetric")
        fdtd.set("simulation time", simulation_time)  
        fdtd.set("use early shutoff", 0)
        fdtd.set("mesh refinement", "conformal variant 1")


        #adding the mode source 
        fdtd.addmode() 
        fdtd.set("name", "TE_Mode") 
        fdtd.set("mode selection", "fundamental mode")  
        fdtd.set("x", source_x) 
        fdtd.set("y min", source_y_min) 
        fdtd.set("y max", source_y_max) 
        fdtd.set("z min", source_z_min) 
        fdtd.set("z max", source_z_max) 
        fdtd.set("center wavelength" , wvl) 
        fdtd.set("wavelength span", bandwidth) 
        # addding the Reflectance monitor
        fdtd.addpower() 
        fdtd.set("name", "Reflectance") 
        fdtd.set("monitor type", "2D X-normal") 
        fdtd.set("x", R_monitor_x) 
        fdtd.set("y min", R_monitor_y_min) 
        fdtd.set("y max", R_monitor_y_max) 
        fdtd.set("z min", R_monitor_z_min) 
        fdtd.set("z max", R_monitor_z_max) 
        fdtd.set("override global monitor settings", 1) 
        fdtd.set("frequency points", freq_res) 
        #adding the Transmittance monitor
        fdtd.addpower() 
        fdtd.set("name", "Transmittance") 
        fdtd.set("monitor type", "2D X-normal") 
        fdtd.set("x",       fdtd_x_max-1e-6) 
        print("T monitor pos=", fdtd_x_max-1e-6)
        fdtd.set("y min", T_monitor_y_min) 
        fdtd.set("y max", T_monitor_y_max) 
        fdtd.set("z min", T_monitor_z_min) 
        fdtd.set("z max", T_monitor_z_max) 
        fdtd.set("override global monitor settings", 1) 
        fdtd.set("frequency points", freq_res) 

        return
# export_simulation_views
def plot_simulation_layout(stru_type, merged_polygon, H, pad_rect_len,bare_mirror_len, source_x, R_monitor_x,
                           T_monitor_x, phc_type, mod_rate, cavity_rect_len=0, wg_thick=0):

    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    # ---- Common conversions ----
    um = 1e6  # meters to microns

    # ---- TOP PLOT: XY view (top-down) ----
    ax = axes[0]

    if stru_type == "no_cav":
    # Left pad
        ax.add_patch(patches.Rectangle((-pad_rect_len * um, -H/2 * um),pad_rect_len * um, H * um,facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5, label='Pad Rect'))

        # Left mirror polygon
        poly_x = merged_polygon[:, 0] * um
        poly_y = merged_polygon[:, 1] * um
        ax.fill(poly_x, poly_y, color='steelblue', alpha=0.6, label='PhC Mirror')
        ax.plot(poly_x, poly_y, color='navy', linewidth=0.5)

        # Right pad
        ax.add_patch(patches.Rectangle(
            (bare_mirror_len * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5))

    else:
        # ---- CAVITY STRUCTURE ----
        # 1. Left mirror left pad
        ax.add_patch(patches.Rectangle(
            (-pad_rect_len * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5, label='Pad Rect'))

        # 2. Left mirror polygon
        poly_x = merged_polygon[:, 0] * um
        poly_y = merged_polygon[:, 1] * um
        ax.fill(poly_x, poly_y, color='steelblue', alpha=0.6, label='Left Mirror')
        ax.plot(poly_x, poly_y, color='navy', linewidth=0.5)

        # 3. Left mirror right pad
        ax.add_patch(patches.Rectangle(
            (bare_mirror_len * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5))

        # 4. Extra pad rect
        extra_pad_start = total_cells * a + pad_rect_len
        ax.add_patch(patches.Rectangle(
            (extra_pad_start * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5))

        # 5. Cavity rect
        end_of_left_mirror = total_cells * a + 2 * pad_rect_len
        cavity_end = end_of_left_mirror + cavity_rect_len
        ax.add_patch(patches.Rectangle(
            (end_of_left_mirror * um, -H/2 * um),
            cavity_rect_len * um, H * um,
            facecolor='gold', alpha=0.5,
            edgecolor='darkorange', linewidth=1, label='Cavity'))

        # 6. Right mirror left pad
        ax.add_patch(patches.Rectangle(
            (cavity_end * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5))

        # 7. Right mirror polygon (flipped)
        right_polygon = merged_polygon.copy()
        right_polygon[:, 0] = bare_mirror_len - right_polygon[:, 0]
        right_poly_x = right_polygon[:, 0] * um + (cavity_end + pad_rect_len) * um
        right_poly_y = right_polygon[:, 1] * um
        ax.fill(right_poly_x, right_poly_y, color='salmon', alpha=0.6, label='Right Mirror')
        ax.plot(right_poly_x, right_poly_y, color='darkred', linewidth=0.5)

        # 8. Right mirror right pad
        end_of_structure = cavity_end + pad_rect_len + total_cells * a
        ax.add_patch(patches.Rectangle(
            (end_of_structure * um, -H/2 * um),
            pad_rect_len * um, H * um,
            facecolor='cornflowerblue', alpha=0.4,
            edgecolor='navy', linewidth=0.5))

        # ---- Source and monitor lines ----
        ax.axvline(x=source_x * um, color='red', linestyle='--',
        linewidth=1.5, label='Mode Source')
        ax.axvline(x=R_monitor_x * um, color='green', linestyle='-.',
        linewidth=1.5, label='R Monitor')
        ax.axvline(x=T_monitor_x * um, color='orange', linestyle='-.',
        linewidth=1.5, label='T Monitor')

        # Labels on lines
        y_label = H/2 * um * 2.5
        ax.text(source_x * um, y_label, 'Source', color='red',
        fontsize=8, ha='center', rotation=90)
        ax.text(R_monitor_x * um, y_label, 'R Mon', color='green',
        fontsize=8, ha='center', rotation=90)
        ax.text(T_monitor_x * um, y_label, 'T Mon', color='orange',
        fontsize=8, ha='center', rotation=90)

        ax.set_xlabel('x (μm)', fontsize=11)
        ax.set_ylabel('y (μm)', fontsize=11)
        ax.set_title(f'XY View — {phc_type} | H={H*1e9:.0f}nm | A={mod_rate*1e9:.0f}nm', fontsize=12)
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=7, ncol=2)
        ax.grid(True, alpha=0.2)

    # ---- BOTTOM PLOT: XZ view (side view) ----
    ax2 = axes[1]

    if stru_type == "no_cav":
    # Full structure as one rect in XZ
        x_start = -pad_rect_len * um
        x_end = (bare_mirror_len + pad_rect_len) * um
        ax2.add_patch(patches.Rectangle(
            (x_start, -wg_thick/2 * um),
            x_end - x_start, wg_thick * um,
            facecolor='steelblue', alpha=0.5,
            edgecolor='navy', linewidth=0.5, label='Waveguide'))
    else:
        x_start = -pad_rect_len * um
        x_end = (end_of_structure + pad_rect_len) * um
        ax2.add_patch(patches.Rectangle(
            (x_start, -wg_thick/2 * um),
            x_end - x_start, wg_thick * um,
            facecolor='steelblue', alpha=0.5,
            edgecolor='navy', linewidth=0.5, label='Waveguide'))

        # Highlight cavity in XZ
        ax2.add_patch(patches.Rectangle(
        (end_of_left_mirror * um, -wg_thick/2 * um),
        cavity_rect_len * um, wg_thick * um,
        facecolor='gold', alpha=0.5,
        edgecolor='darkorange', linewidth=1, label='Cavity'))

        # Same source/monitor lines
        ax2.axvline(x=source_x * um, color='red', linestyle='--', linewidth=1.5)
        ax2.axvline(x=R_monitor_x * um, color='green', linestyle='-.', linewidth=1.5)
        ax2.axvline(x=T_monitor_x * um, color='orange', linestyle='-.', linewidth=1.5)

        ax2.set_xlabel('x (μm)', fontsize=11)
        ax2.set_ylabel('z (μm)', fontsize=11)
        ax2.set_title('XZ View (Side)', fontsize=12)
        ax2.set_aspect('equal')
        ax2.legend(loc='upper right', fontsize=7)
        ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    fname = f"{phc_type}_{stru_type}_H_{H*1e9:.0f}nm_A_{mod_rate*1e9:.0f}nm_layout.png"
    plt.savefig(fname, dpi=200)
    plt.close()
    print(f"  Layout saved: {fname}")
    return

def loop_simulation(parameter, stru_type, phc_type):
    all_T = []
    all_R = []          # list of reflectance arrays
    all_f = None        # frequency array (same for all runs)
    all_labels = []     # labels for each run

    os.makedirs("plots", exist_ok=True)
    os.makedirs("data", exist_ok=True)

    for i in range(len(parameter)):
        H = H_values[i]
        mod_rate = Amp_values
        bare_mirror_len = total_cells * a
        wg_thick = Z_thk

        print(f"\n{'='*60}")
        print(f"Run {i+1}/{len(parameter)}: H = {H*1e9:.0f} nm, A = {mod_rate*1e9:.0f} nm")
        print(f"{'='*60}")

        # Build geometry
        make_fdtd_simu_domain(stru_type,phc_type, H, mod_rate, wg_thick, pad_rect_len, bare_mirror_len)

        # Plot layout check
        merged_polygon = make_phc_(phc_type, mod_rate,H)
        T_monitor_x = fdtd_x_max
        
        plot_simulation_layout(
            stru_type, merged_polygon, H, pad_rect_len,
            bare_mirror_len, source_x, R_monitor_x, T_monitor_x,
            phc_type, mod_rate,
            cavity_rect_len=cavity_rect_len,
            wg_thick=Z_thk
        )

        # Save and run
        run_name = f"{phc_type}_{stru_type}_H_{H*1e9:.0f}nm_A_{mod_rate*1e9:.0f}nm"
        fdtd.save(f"data/{run_name}.fsp")
        fdtd.run()

        # Get data
        T = np.squeeze(fdtd.transmission("Transmittance"))
        R = np.squeeze(fdtd.transmission("Reflectance"))
        f = np.squeeze(fdtd.getdata("Transmittance", "f"))
        wavelength = 3e8 / f * 1e9   # convert to nm

        # Store
        all_T.append(T)
        all_R.append(R)
        if all_f is None:
            all_f = f
        all_labels.append(f"H={H*1e9:.0f}nm")

        # ---- Plot individual wavelength vs T ----
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(wavelength, T, 'b-', linewidth=1)
        ax.set_xlabel('Wavelength (nm)', fontsize=12)
        ax.set_ylabel('Transmission', fontsize=12)
        ax.set_title(f'{phc_type} | H = {H*1e9:.0f} nm | A = {mod_rate*1e9:.0f} nm', fontsize=13)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"plots/{run_name}_transmission.png", dpi=200)
        plt.close()
        print(f"  Plot saved: plots/{run_name}_transmission.png")

        # Save individual run data via matlabsave
        fdtd.matlabsave(os.path.abspath(f"data/{run_name}"))
        print(f"  Data saved: data/{run_name}.mat")

        # ---- Switch to layout mode and clear for next run ----
        if i < len(parameter) - 1:
            fdtd.switchtolayout()
            fdtd.deleteall()
            print("  Cleared layout for next run")

    # ============================================
    # After all runs: combined plot
    # ============================================
    wavelength = 3e8 / all_f * 1e9

    fig, ax = plt.subplots(figsize=(12, 6))
    for idx in range(len(all_T)):
        ax.plot(wavelength, all_T[idx], linewidth=1, label=all_labels[idx])

    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Transmission', fontsize=12)
    ax.set_title(f'{phc_type} — All Cavity Resonances', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"plots/{phc_type}_{stru_type}_all_resonances.png", dpi=200)
    plt.close()
    print(f"\nCombined plot saved: plots/{phc_type}_{stru_type}_all_resonances.png")

    # ============================================
    # Save all cavity resonances as one .mat file
    # ============================================
    # Put the arrays into Lumerical workspace and save
    fdtd.putv("all_T", np.array(all_T))
    fdtd.putv("all_R", np.array(all_R))
    fdtd.putv("frequency", all_f)
    fdtd.putv("wavelength_nm", wavelength)
    fdtd.matlabsave(f"data/{phc_type}_{stru_type}_all_cavity_resonances")
    print(f"All data saved: data/{phc_type}_{stru_type}_all_cavity_resonances.mat")

    return all_T, all_R, all_f

# Run it
all_T, all_R, all_f = loop_simulation(H_values, "cavity", "fence_gate")
input("Press Enter to close the session")









