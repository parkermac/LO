"""
Evaluation plotting script for Penn Cove mussel raft sonde data.

Produces two figure types per station per year:
  1. Heatmap (time × depth) — overview of the full record.
  2. Time series (one line per depth) — fine-grained inspection per variable.

Usage: set `year` at the top and run.  Set year = None to loop all years.

Finalized for public use: 2026/04/13

Author: Dakota Mascarenas
"""

# %%

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from lo_tools import Lfun

Ldir = Lfun.Lstart()

# ---- USER SETTINGS -------------------------------------------------------
year = None   # int (e.g. 2016) or None to loop through all available years

# Optional y-axis zoom per variable (overrides vmin/vmax for display only).
# Set to None to use the default range from var_meta.
# Example: {'CT': (8, 14), 'SA': (28, 32)}
ylim_override = {
    'CT':           None,
    'SA':           None,
    'DO (uM)':      (0, 600),
    'PH':           (7.2, 8.8),
    'Chl (mg m-3)': (0, 300),
}
# --------------------------------------------------------------------------

source = 'pcRaft'
otype  = 'moor'
in_dir = Ldir['LOo'] / 'obs' / source / otype

plot_dir = Ldir['LOo'] / 'obs' / source / 'eval_plots'
Lfun.make_dir(plot_dir)

# Variable display metadata: name -> (colormap, (vmin, vmax), label)
var_meta = {
    'CT':          ('RdYlBu_r', (0, 25),    'Cons. Temp. (°C)'),
    'SA':          ('viridis',  (0, 33),   'Abs. Salinity (g kg⁻¹)'),
    'DO (uM)':     ('YlGnBu',  (0, 600),   'DO (μM)'),
    'PH':          ('PiYG',     (7.2, 8.8), 'pH'),
    'Chl (mg m-3)':('YlGn',    (0, 300),    'Chl (mg m⁻³)'),
}

stations = ['raft_main', 'raft_hiRes']

# %%

def plot_year(ds, station, yr, plot_dir):
    """
    Create a heatmap figure (time × depth) for each variable in ds,
    restricted to the given year.  Saves one PNG per station per year.
    """

    # Slice to the requested year
    time_idx = pd.DatetimeIndex(ds['time'].values)
    mask = time_idx.year == yr
    if mask.sum() == 0:
        print(f'  No data for {station} in {yr}, skipping.')
        return

    ds_yr = ds.isel(time=np.where(mask)[0])
    time_vals = pd.DatetimeIndex(ds_yr['time'].values)
    z_vals    = ds_yr['z'].values   # negative depths (e.g. -7 to -0.5)

    avail_vars = [v for v in var_meta if v in ds_yr.data_vars]
    n_vars = len(avail_vars)
    if n_vars == 0:
        print(f'  No plottable variables in {station} for {yr}.')
        return

    fig, axes = plt.subplots(n_vars, 1,
                             figsize=(14, 3.5 * n_vars),
                             sharex=True)
    if n_vars == 1:
        axes = [axes]

    fig.suptitle(f'{station} — {yr}', fontsize=13, fontweight='bold')

    for ax, vname in zip(axes, avail_vars):
        cmap, (vmin, vmax), label = var_meta[vname]
        # Apply display zoom if set
        override = ylim_override.get(vname)
        disp_min, disp_max = override if override is not None else (vmin, vmax)
        data = ds_yr[vname].values   # shape (NT, NZ)

        # pcolormesh needs grid edges; build them from centres
        # time axis: use numeric values (matplotlib dates)
        t_num  = mdates.date2num(time_vals.to_pydatetime())
        dt     = np.diff(t_num)
        dt     = np.append(dt, dt[-1])
        t_edges = np.append(t_num - dt / 2, t_num[-1] + dt[-1] / 2)

        dz      = np.diff(z_vals)
        dz      = np.append(dz[0], dz)
        z_edges = np.append(z_vals - np.abs(dz) / 2,
                             z_vals[-1] + np.abs(dz[-1]) / 2)

        # Use masked array so NaN shows as white
        data_masked = np.ma.masked_invalid(data.T)   # (NZ, NT)

        pcm = ax.pcolormesh(t_edges, z_edges, data_masked,
                    cmap=cmap, vmin=disp_min, vmax=disp_max,
                    shading='flat')

        cb = fig.colorbar(pcm, ax=ax, pad=0.01, fraction=0.02)
        cb.set_label(label, fontsize=8)

        # Scatter any values outside the expected range so they stand out
        out_mask = (~np.isnan(data)) & ((data < vmin) | (data > vmax))
        if out_mask.any():
            t_idx, z_idx = np.where(out_mask)
            ax.scatter(t_num[t_idx], z_vals[z_idx],
                       color='red', s=10, zorder=5,
                       label=f'out-of-range ({out_mask.sum()})')
            ax.legend(loc='upper right', fontsize=7)

        ax.set_ylabel('Depth (m)', fontsize=9)
        ax.set_title(label, fontsize=9)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.grid(True, alpha=0.25, linestyle=':')

    axes[-1].set_xlabel(f'{yr} (UTC)', fontsize=9)
    plt.tight_layout()

    out_fn = plot_dir / f'{station}_{yr}_heatmap.png'
    plt.savefig(out_fn, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_fn}')


def plot_year_timeseries(ds, station, yr, plot_dir):
    """
    Create a time-series figure for each variable in ds, restricted to the
    given year.  Each depth level is drawn as a separate line.  Values outside
    the expected range are highlighted with red markers.
    Saves one PNG per station per year.
    """

    time_idx = pd.DatetimeIndex(ds['time'].values)
    mask = time_idx.year == yr
    if mask.sum() == 0:
        return  # already reported by plot_year

    ds_yr = ds.isel(time=np.where(mask)[0])
    time_vals = pd.DatetimeIndex(ds_yr['time'].values)
    z_vals    = ds_yr['z'].values

    avail_vars = [v for v in var_meta if v in ds_yr.data_vars]
    n_vars = len(avail_vars)
    if n_vars == 0:
        return

    # Use a colormap to assign one colour per depth level
    depth_cmap = plt.get_cmap('plasma', max(len(z_vals), 1))

    fig, axes = plt.subplots(n_vars, 1,
                             figsize=(14, 3.5 * n_vars),
                             sharex=True)
    if n_vars == 1:
        axes = [axes]

    fig.suptitle(f'{station} — {yr} (time series)', fontsize=13, fontweight='bold')

    for ax, vname in zip(axes, avail_vars):
        _, (vmin, vmax), label = var_meta[vname]
        # Apply display zoom if set
        override = ylim_override.get(vname)
        disp_min, disp_max = override if override is not None else (vmin, vmax)
        data = ds_yr[vname].values   # shape (NT, NZ)

        for z_idx, z in enumerate(z_vals):
            col = depth_cmap(z_idx / max(len(z_vals) - 1, 1))
            series = data[:, z_idx]
            valid  = ~np.isnan(series)
            if valid.sum() == 0:
                continue
            ax.plot(time_vals[valid], series[valid],
                    color=col, linewidth=0.8, alpha=0.85,
                    label=f'z={z:.1f} m')

            # Highlight out-of-range points
            out = valid & ((series < vmin) | (series > vmax))
            if out.any():
                ax.scatter(time_vals[out], series[out],
                           color='red', s=14, zorder=5)

        # Shaded band for expected range
        ax.axhspan(vmin, vmax, color='lightgrey', alpha=0.15, zorder=0)
        ax.set_ylim(disp_min, disp_max)
        ax.set_ylabel(label, fontsize=9)
        ax.set_title(label, fontsize=9)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.grid(True, alpha=0.25, linestyle=':')
        ax.legend(loc='upper right', fontsize=7, ncol=2)

    axes[-1].set_xlabel(f'{yr} (UTC)', fontsize=9)
    plt.tight_layout()

    out_fn = plot_dir / f'{station}_{yr}_timeseries.png'
    plt.savefig(out_fn, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_fn}')


# %%

for sn in stations:
    nc_path = in_dir / f'{sn}.nc'
    if not nc_path.exists():
        print(f'File not found: {nc_path}')
        continue

    ds = xr.open_dataset(nc_path)

    all_years = sorted(set(pd.DatetimeIndex(ds['time'].values).year))
    years_to_plot = [year] if year is not None else all_years

    print(f'\n{sn}: available years = {all_years}')
    for yr in years_to_plot:
        print(f'  Plotting {yr}...')
        plot_year(ds, sn, yr, plot_dir)
        plot_year_timeseries(ds, sn, yr, plot_dir)

    ds.close()

print('\nDone.')
