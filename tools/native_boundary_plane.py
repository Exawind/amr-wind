# -*- coding: utf-8 -*-

"""\
AMR-Wind Native Boundary Plane
------------------------------

This module provides a python class :class:`NativeBoundaryPlane` that can be used to define, manipulate, and output native boundary planes for AMR-Wind

"""
from amrex_plotfile import AmrexPlotFile
import amrex.space3d as amr
import amrex_utils as au


class NativeBoundaryPlane:
    def __init__(self, iname, field, ncomp, ori, step, time, odir):

        self.odir = odir
        hname = self.odir / f"Header_{ori}_{field}"

        spacedim = amr.Config.spacedim
        pp = amr.ParmParse("")
        pp.addfile(iname)
        pp_geom = amr.ParmParse("geometry")
        prob_lo = [
            pp_geom.get_real("prob_lo", 0),
            pp_geom.get_real("prob_lo", 1),
            pp_geom.get_real("prob_lo", 2),
        ]
        prob_hi = [
            pp_geom.get_real("prob_hi", 0),
            pp_geom.get_real("prob_hi", 1),
            pp_geom.get_real("prob_hi", 2),
        ]
        pp_amr = amr.ParmParse("amr")
        n_cell = [
            pp_amr.get_int("n_cell", 0),
            pp_amr.get_int("n_cell", 1),
            pp_amr.get_int("n_cell", 2),
        ]
        finest_level = pp_amr.get_int("max_level")
        nlevels = finest_level + 1
        ref_ratio = [2] * (nlevels - 1)
        vnames = au.variable_names(field, ncomp)

        normal = au.normal_from_ori(ori)
        perp = au.perpendicular_from_ori(ori)
        domain_lo = amr.IntVect(0)
        domain_hi = amr.IntVect(n_cell[0] - 1, n_cell[1] - 1, n_cell[2] - 1)
        domain_lo[normal] = -1 if au.is_low(ori) else domain_hi[normal]
        domain_hi[normal] = 0 if au.is_low(ori) else domain_hi[normal] + 1
        domain = amr.Box(domain_lo, domain_hi)

        dx_normal = (prob_hi[normal] - prob_lo[normal]) / (n_cell[normal])
        plo_normal = prob_lo[normal]
        phi_normal = prob_hi[normal]
        n_cell[normal] = 2
        prob_lo[normal] = (
            plo_normal - dx_normal if au.is_low(ori) else phi_normal - dx_normal
        )
        prob_hi[normal] = (
            plo_normal + dx_normal if au.is_low(ori) else phi_normal + dx_normal
        )

        cell_sizes = []
        for ilev in range(nlevels):
            cs = [(prob_hi[i] - prob_lo[i]) / (n_cell[i]) for i in range(spacedim)]
            if ilev > 0:
                cs = [x / (ref_ratio[ilev - 1] ** ilev) for x in cs]
            cell_sizes.append(cs)

        prob_domain = [domain]
        bas = [amr.BoxArray(domain)]
        for i in range(1, nlevels):
            small_end = prob_domain[ilev - 1].small_end
            big_end = prob_domain[ilev - 1].big_end
            for p in perp:
                small_end[p] *= ref_ratio[ilev - 1]
                big_end[p] = (big_end[p] + 1) * ref_ratio[ilev - 1] - 1
            if not au.is_low(ori):
                small_end[normal] = small_end[normal] * ref_ratio[ilev - 1] + 1
                big_end[normal] = small_end[normal] + 1

            prob_domain.append(amr.Box(small_end, big_end))
            bas.append(amr.BoxArray(prob_domain[ilev]))

        mfs = []
        n_grow = 0
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

        self.plt = AmrexPlotFile(hname)
        lvl_pfx = "Level_"
        mf_names = [f"{lvl_pfx}{x}/{field}_{ori}" for x in range(nlevels)]
        self.plt.define(
            mfs,
            vnames,
            mf_names,
            step,
            time,
            prob_lo,
            prob_hi,
            ref_ratio,
            prob_domain,
            glohis,
        )

    def evaluate(self, functors):
        self.plt.evaluate(functors)

    def write(self):
        overwrite = True
        self.plt.write(self.odir, overwrite)
