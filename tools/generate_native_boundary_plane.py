"""\
A tool to generate native boundary planes
-----------------------------------------

This script is for generating boundary plane data.

"""
import pathlib
import numpy as np
import amrex.space3d as amr
from amrex_plotfile import AmrexPlotFile
import amrex_utils as au

def main():
    amr.initialize([])

    bdir = pathlib.Path("bndry_file")
    
    field = "velocity"
    ori = 0
    ncomp = 3
    spacedim = amr.Config.spacedim
    lvl_pfx = "Level_"

    # fixme parm parse these
    prob_lo = [0, 0, 0]
    prob_hi = [1000, 1000, 1000]
    n_cell = [48, 48, 48]
    nlevels = 2
    ref_ratio = [2]
    n_grow = 0

    vnames = au.variable_names(field, ncomp)

    # adjust for plane
    normal = au.normal_from_ori(ori)
    perp = au.perpendicular_from_ori(ori)
    domain_lo = amr.IntVect(0)
    domain_hi = amr.IntVect(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)
    domain_lo[normal] = -1 if au.is_low(ori) else domain_hi[normal]
    domain_hi[normal] = 0 if au.is_low(ori) else domain_hi[normal]+1
    domain = amr.Box(domain_lo, domain_hi)

    dx_normal = (prob_hi[normal] - prob_lo[normal]) / (n_cell[normal])
    plo_normal = prob_lo[normal]
    phi_normal = prob_hi[normal]
    n_cell[normal] = 2
    prob_lo[normal] = plo_normal - dx_normal if au.is_low(ori) else phi_normal - dx_normal
    prob_hi[normal] = plo_normal + dx_normal if au.is_low(ori) else phi_normal + dx_normal

    cell_sizes = []
    for ilev in range(nlevels):
        cs = [(prob_hi[i] - prob_lo[i]) / (n_cell[i]) for i in range(spacedim)]
        if ilev > 0:
            cs = [x/(ref_ratio[ilev-1] ** ilev) for x in cs]
        cell_sizes.append(cs)

    prob_domain = [domain]
    bas = [amr.BoxArray(domain)]
    for i in range(1, nlevels):
        small_end = prob_domain[ilev-1].small_end
        big_end = prob_domain[ilev-1].big_end
        for p in perp:
            small_end[p] *= ref_ratio[ilev-1]
            big_end[p] = (big_end[p] + 1) * ref_ratio[ilev-1] - 1
        if not au.is_low(ori):
            small_end[normal] = small_end[normal] * ref_ratio[ilev-1] + 1
            big_end[normal] = small_end[normal]+1
            
        prob_domain.append(amr.Box(small_end, big_end))
        bas.append(amr.BoxArray(prob_domain[ilev]))
    
    # time/step data
    nsteps = 100
    steps = [x for x in range(nsteps+1)]
    times = np.linspace(0, 5, nsteps+1)

    for step, time in zip(steps, times):
        odir = bdir / f"bndry_output{step:05d}"
        hname = odir / f"Header_{ori}_{field}"
        print(step, time, hname)

        mfs = []
        for ilev in range(nlevels):
            dm = amr.DistributionMapping(bas[ilev])
            mf = amr.MultiFab(bas[ilev], dm, ncomp, n_grow)
            mf.set_val(0.0)
            mfs.append(mf)

        glohis = [[] for _ in range(nlevels)]
        for ilev in range(nlevels):
            nglohis = []
            for igrid in range(mfs[ilev].box_array().size):
                glos = []
                ghis = []
                for idim in range(spacedim):
                    if idim == normal:
                        offset = (prob_lo[idim] + prob_hi[idim]) / 2
                        glo = offset - cell_sizes[ilev][idim]
                        ghi = offset + cell_sizes[ilev][idim]
                    else:
                        glo = prob_lo[idim]
                        ghi = prob_hi[idim]
                    glos.append(glo)
                    ghis.append(ghi)
                nglohis.append([glos, ghis])
            glohis[ilev] = nglohis

        mf_names = [f"{lvl_pfx}{x}/{field}_{ori}" for x in range(nlevels)]

        plt = AmrexPlotFile(hname)
        plt.define(mfs, vnames, mf_names, step, time, prob_lo, prob_hi, ref_ratio, prob_domain, glohis)

        for ilev in range(nlevels):
            data = plt.mfs[ilev].to_xp()[0]
            los = plt.glohis[ilev][0][0]
            his = plt.glohis[ilev][0][1]
            xlo = los[perp[0]]
            xhi = his[perp[0]]
            ylo = los[perp[1]]
            yhi = his[perp[1]]
            nx = plt.mfs[ilev].box_array()[0].size[perp[0]]
            ny = plt.mfs[ilev].box_array()[0].size[perp[1]]
            dx = plt.cell_sizes[ilev][perp[0]]
            dy = plt.cell_sizes[ilev][perp[1]]
            x = np.linspace(xlo + 0.5 * dx, xhi - 0.5 * dx, nx)
            y = np.linspace(ylo + 0.5 * dy, yhi - 0.5 * dy, ny)
            xg, yg = np.meshgrid(x, y, indexing="ij")
            xc = 500
            sigma = 100
            for i in range(data.shape[normal]):
                for nc in range(plt.ncomp):
                    slc = au.slice_from_normal(normal, i, nc)
                    data[slc] = np.exp(-((xg-xc)**2/(2*sigma**2) + (yg-xc)**2/(2*sigma**2))) * time

        plt.write(hname.parent)

    np.savetxt(bdir / "time.dat", np.c_[steps, times], fmt="%.17g")
    
if __name__ == "__main__":
    main()

