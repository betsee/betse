#!/usr/bin/env python3
# ....................{ LICENSE                           }....................
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................

import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import colorbar
from matplotlib import rcParams
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Circle
from matplotlib import path

from betse.util.math.geometry.polygon.geopolyconvex import clip_counterclockwise
from betse.util.math.geometry.polygon.geopoly import orient_counterclockwise
from betse.util.io.log import logs

from scipy.spatial import cKDTree, Delaunay

# FIXME get face to edge mappings for tri and vor!
# FIXME make mids mappers
class DECMesh(object):
    """
    Creates a Discrete Exterior Calculus mesh system, with primal triangulation
    and dual Voronoi meshes from a set of seed points.

    """

    def __init__(self,
                 cell_radius,
                 seed_points=None,
                 alpha_shape=0.4,
                 use_alpha_shape=True,
                 single_cell_noise=0.5,
                 single_cell_sides=6,
                 image_mask=None,
                 make_all_operators = True):

        self.single_cell_noise = single_cell_noise
        self.single_cell_sides = single_cell_sides
        self.alpha_shape = alpha_shape
        self.cell_radius = cell_radius
        self.use_alpha_shape = use_alpha_shape
        self.image_mask = image_mask
        self.make_all_operators = make_all_operators

        if seed_points is None:
            self.make_single_cell_points()
        else:
            self.tri_verts = seed_points


        self.init_mesh()

    def init_mesh(self):

        self.create_tri_mesh()
        self.create_mappings()
        self.process_vormesh()
        self.create_core_operators()

        if self.make_all_operators:
            self.create_aux_operators()

    def make_single_cell_points(self):
        """
        Creates seed points (tri_mesh verts) for a single closed voronoi cell and
        its boundary.

        """

        centre = [0.0, 0.0]

        angles = [(2 * n * np.pi) / self.single_cell_sides
                  for n in range(self.single_cell_sides)]

        xenv = 2 * self.cell_radius * np.cos(angles)
        yenv = 2 * self.cell_radius * np.sin(angles)

        noisex = np.random.random(self.single_cell_sides + 1)
        noisey = np.random.random(self.single_cell_sides + 1)

        xenv = np.hstack((xenv, 0.0)) + noisex * self.cell_radius * self.single_cell_noise
        yenv = np.hstack((yenv, 0.0)) + noisey * self.cell_radius * self.single_cell_noise

        self.tri_verts = np.column_stack((xenv, yenv))

    def create_tri_mesh(self):
        """
        Calculates the Delaunay triangulation and non-convex Hull (using alpha shapes).

        """
        # calculate the Delaunday triangulation based on the cluster-masked seed points:
        trimesh = Delaunay(self.tri_verts)  # Delaunay trianulation of cell centres


        self.n_tverts = len(self.tri_verts)

        # Calculate the centroids of the triangulation:
        self.tri_circcents = []
        self.vor_verts_bound = [] # vertices at the boundary
        self.tri_rcircs = [] # circumradius of triangle
        self.tri_rin = [] # inradius of triangle
        self.tricell_i = []  # index of simplexes
        self.tri_sa = []  # surface area of triangle faces
        self.tcell_verts = []  # x,y components of tri_mesh cells
        simplices2 = []

        for i, vert_inds in enumerate(trimesh.simplices):

            abc = trimesh.points[vert_inds]
            vx, vy, r_circ, r_in = self.circumc(abc[0], abc[1], abc[2])

            sa = self.area(abc)  # surface area of triangle

            if self.use_alpha_shape:

                # exclude large circumcircles from the triangulation according to "alpha-shape"
                # methodology:
                if r_circ < (self.cell_radius) / self.alpha_shape:

                    # but also check that circumcentre is in the cluster mask...
                    if self.image_mask is not None:
                        flagc = self.image_mask.clipping_function(vx, vy)

                    else:
                        flagc = 1.0

                    if flagc != 0.0:  # if it's not outside the cluster region, include the simplex:
                        self.tri_circcents.append([vx, vy])
                        self.tri_rcircs.append(r_circ)
                        self.tri_rin.append(r_in)
                        self.tricell_i.append(i)
                        simplices2.append(vert_inds)
                        self.tri_sa.append(sa)
                        self.tcell_verts.append(abc)

                self.tri_cells = np.asarray(simplices2)  # reassign point inds to retained simplices

            else:
                self.tri_cells = trimesh.simplices
                self.tri_circcents.append([vx, vy])
                self.tri_rcircs.append(r_circ)
                self.tri_rin.append(r_in)
                self.tricell_i.append(i)
                self.tri_sa.append(sa)
                self.tcell_verts.append(abc)

        self.n_tcell = len(self.tricell_i)  # number of simplexes in trimesh

        # Calculate edges for cells trimesh, and the hull:
        all_edges = set()  # set of all vertex pairs
        unique_edges = set()  # set of unique vertex pairs
        hull_points = []
        hull_edges = []

        for vi, vj, vk in self.tri_cells:
            all_edges.add((vi, vj))
            all_edges.add((vj, vk))
            all_edges.add((vk, vi))

        for va, vb in all_edges:

            if (va, vb) in all_edges and (vb, va) not in all_edges:
                # if there isn't a double-pair, then add these edges to the hull:
                # (this is based on the logic that when traversing the points of the
                # triangular simplices, only the boundary edges are traversed once,
                # since they don't have a neighbouring simplex at the bounds.)
                hull_points.append(va)
                hull_points.append(vb)

                hull_edges.append([va, vb])

                # calculate the boundary vor verts from hull edge midpoint:
                pta = self.tri_verts[va]
                pself = self.tri_verts[vb]

                vorb = (pta + pself) / 2

                self.vor_verts_bound.append(vorb)

            if (vb, va) not in unique_edges:
                unique_edges.add((va, vb))

        self.tri_circcents = np.asarray(self.tri_circcents)
        self.vor_verts_bound = np.asarray(self.vor_verts_bound)
        self.tri_rcircs = np.asarray(self.tri_rcircs)
        self.tri_rin = np.asarray(self.tri_rin)
        self.tri_sa = np.asarray(self.tri_sa)
        self.tri_cell_i = np.asarray(self.tricell_i)
        self.tcell_verts = np.asarray(self.tcell_verts)

        self.tri_vert_i = np.linspace(0, self.n_tverts - 1,
                                      self.n_tverts, dtype=np.int)

        self.bflags_tverts = np.unique(hull_points)
        self.tri_edges = np.asarray(list(unique_edges))
        self.n_tedges = len(self.tri_edges)  # number of edges in trimesh


        # Process edges to create flags of edge indices:
        tri_edges = self.tri_edges.tolist()

        bflags_tedges = []

        for vi, vj in hull_edges:
            if [vi, vj] in tri_edges:
                kk = tri_edges.index([vi, vj])
            elif [vj, vi] in tri_edges:
                kk = tri_edges.index([vj, vi])
            bflags_tedges.append(kk)

        self.bflags_tedges = np.asarray(bflags_tedges) # indices of edges on the boundary

        # inds to full voronoi cells only (official cells of cluster):
        self.biocell_i = np.delete(self.tri_vert_i, self.bflags_tverts)

        self.tri_edge_i = np.linspace(0, self.n_tedges - 1, self.n_tedges, dtype=np.int)

        # include hull-edge-mids in voronoi vert array:
        # self.biocell_verts = self.tri_circcents*1 # assign core vor_verts to new data_structure

        # reassign vor_verts to contain boundary verts:
        self.vor_verts = np.vstack((self.tri_circcents, self.vor_verts_bound))

    def create_mappings(self):


        # make the face-to-edge indices mapping for the tri_mesh:
        tri_edges = self.tri_edges.tolist()

        face_to_edges = [[] for ii in range(self.n_tcell)]
        bflags_tcells = []

        for ci, (vi, vj, vk) in enumerate(self.tri_cells):

            if [vi, vj] not in tri_edges:
                # get the index of the opposite sign edge:
                ea = tri_edges.index([vj, vi])
                face_to_edges[ci].append(ea) # append the edge index to the array at the cell index

            else:
                # get the forward sign edge:
                ea = tri_edges.index([vi, vj])
                face_to_edges[ci].append(ea)

            if [vj, vk] not in tri_edges:
                # get the index of the opposite sign edge:
                eb = tri_edges.index([vk, vj])
                face_to_edges[ci].append(eb)
            else:
                # get the forward sign edge:
                eb = tri_edges.index([vj, vk])
                face_to_edges[ci].append(eb)

            if [vk, vi] not in tri_edges:
                # get the index of the opposite sign edge:
                ec = tri_edges.index([vi, vk])
                face_to_edges[ci].append(ec)
            else:
                # get the forward sign edge:
                ec = tri_edges.index([vk, vi])
                face_to_edges[ci].append(ec)

            # if any edge is on the boundary, mark the cell
            if ea in self.bflags_tedges or eb in self.bflags_tedges or ec in self.bflags_tedges:
                bflags_tcells.append(ci)

        self.tface_to_tedges = np.asarray(face_to_edges) # tri_face index to tri_edges indices mapping
        self.bflags_tcells = np.asarray(bflags_tcells)

        # create an array giving a list of simplex indices for each tri_vert
        verts_to_simps = [[] for i in range(len(self.tri_verts))]

        for si, vert_inds in enumerate(self.tri_cells):

            for vi in vert_inds:
                verts_to_simps[vi].append(si)

        self.tverts_to_tcell = np.asarray(verts_to_simps)

        # For each trivert, get indices from vorvert to complete Voronoi cells
        # on the boundary:
        extra_vor_inds = [[] for i in self.tri_verts]

        for ii, ei in enumerate(self.bflags_tedges):

            vi, vj = self.tri_edges[ei]
            # get the vorvert corresponding to the hull_edge:
            ptv = self.vor_verts_bound[ii]

            # get the indices to the vorvert in vor_verts:
            indv = np.where(self.vor_verts == ptv)[0][0]

            # append the shared vor vert ind to both points:
            extra_vor_inds[vi].append(indv)
            extra_vor_inds[vj].append(indv)

        extra_vor_inds = np.asarray(extra_vor_inds)

        # create an array giving a list of vorvert indices from vor_verts for each tri_vert
        triverts_to_vorverts = [[] for i in range(len(self.tri_verts))]

        vorpoints_tree = cKDTree(self.vor_verts)

        for vi, (simp_inds, extra_inds) in enumerate(zip(self.tverts_to_tcell,
                                                         extra_vor_inds)):

            vpts = self.vor_verts[simp_inds]
            _, all_vort_inds = vorpoints_tree.query(vpts)

            for i in all_vort_inds:
                triverts_to_vorverts[vi].append(i)

            # if there's anything in the extra inds array, add both partners
            # to the verts defining the Voronoi region:
            if len(extra_inds):
                triverts_to_vorverts[vi].append(extra_inds[0])
                triverts_to_vorverts[vi].append(extra_inds[1])

        self.triverts_to_vorverts = np.asarray(triverts_to_vorverts)

    def process_vormesh(self):

        vor_cells = []  # vor_vert inds specifying dual cell for each tri_vert
        vcell_verts = []
        vor_cents = []  # centroid of voronoi polygons
        vor_sa = []  # surface area of voronoi polygons

        for vi, vort_inds in enumerate(self.triverts_to_vorverts):
            vort_inds = np.asarray(vort_inds)

            # get voronoi verts corresponding to each simplex:
            vor_pts = self.vor_verts[vort_inds]

            cent = vor_pts.mean(axis=0)  # calculate the centre point
            angles = np.arctan2(vor_pts[:, 1] - cent[1], vor_pts[:, 0] - cent[0])  # calculate point angles
            sorted_region = vort_inds[np.argsort(angles)]  # sort indices counter-clockwise
            vor_cells.append(sorted_region)  # add sorted list to the regions structure

            vpts = self.vor_verts[sorted_region]  # x,y coordinates, sorted
            vcell_verts.append(vpts)

            # Calculate centroid and area of the voronoi polygon:
            cx, cy = self.poly_centroid(vpts)
            vor_cents.append([cx, cy])

            sa = self.area(vpts)
            vor_sa.append(sa)

        self.vor_cells = np.asarray(vor_cells)
        self.vcell_verts = np.asarray(vcell_verts)

        self.vor_cents = np.asarray(vor_cents)
        self.vor_sa = np.asarray(vor_sa)

        # find edges of Voronoi dual mesh:
        # Want vor edges to have the same index as tri_edges and to be
        # perpendicular bisectors

        # have one vor_edge vert pair (of all_vor_inds) for each tri-edge
        vor_edges = []

        tri_tang = []  # tangent vectors to tri_edges
        vor_tang = []  # tangent vectors to vor_edges

        vor_edge_len = []  # length of vor_edge
        tri_edge_len = []  # length of tri_edge

        #         vor_edge_mids = [] # mids of Voronoi cell edges (not needed as same as tri!)
        tri_edge_mids = []  # mids of tri cell edges

        for tei, (ti, tj) in enumerate(self.tri_edges):

            # get the voronoi indices (each sorted counter clockwise) for
            # defining voronoi dual cells surrounding each tri vert point:
            vor_reg_i = self.triverts_to_vorverts[ti]
            vor_reg_j = self.triverts_to_vorverts[tj]

            # find the (hopefully 2!) common inds between the two cells,
            # representing the shared edge:
            shared_ij = np.intersect1d(vor_reg_i, vor_reg_j)

            # assert (len(shared_ij) == 2), "Shared vor_cell edge inds must be length 2"

            if len(shared_ij) == 2:

                # find points representing vor and tri edges:

                tpi = self.tri_verts[ti]
                tpj = self.tri_verts[tj]

                tan_t = tpj - tpi

                tri_len = np.linalg.norm(tan_t)
                #     tri_len = np.sqrt(tan_t[0]**2 + tan_t[1]**2)

                pptm = (tpi + tpj) / 2  # calculate midpoint of tri-edge:

                tri_edge_mids.append(pptm)

                assert (tri_len != 0.0), "Tri-edge length equal to zero!"

                tri_edge_len.append(tri_len * 1)
                tri_tang.append(1 * tan_t / tri_len)

                vpi = self.vor_verts[shared_ij[0]]
                vpj = self.vor_verts[shared_ij[1]]

                tan_v = vpj - vpi

                vor_len = np.linalg.norm(tan_v)
                #     vor_len = np.sqrt(tan_v[0]**2 + tan_v[1]**2)

                #             ppvm = (vpi + vpj)/2 # midpoint of voronoi edge
                #             vor_edge_mids.append(ppvm)

                assert (vor_len != 0.0), "Vor-edge length equal to zero!"

                vor_edge_len.append(vor_len * 1)

                cross_tp = tan_t[0] * tan_v[1] - tan_v[0] * tan_t[1]

                # Check that the tri_edge and vor_edge are orthogonal:
                dottv = np.dot(tan_t, tan_v)
                dot_check = np.round(np.abs(dottv), 20)

                # if yes, add vor_edge points with 90 degree clockwise rotation to
                # the tri_edge:
                if dot_check == 0.0:

                    if np.sign(cross_tp) == -1.0:

                        vor_edges.append([shared_ij[0], shared_ij[1]])
                        vor_tang.append(tan_v / vor_len)

                    elif np.sign(cross_tp) == 1.0:

                        vor_edges.append([shared_ij[1], shared_ij[0]])
                        vor_tang.append(-tan_v / vor_len)

        self.vor_edges = np.asarray(vor_edges)
        self.tri_tang = np.asarray(tri_tang)
        self.vor_tang = np.asarray(vor_tang)

        self.vor_edge_len = np.asarray(vor_edge_len)  # length of vor_edge
        self.tri_edge_len = np.asarray(tri_edge_len)  # length of tri_edge

        self.tri_mids = np.asarray(tri_edge_mids)  # edge midpoints (same vals for vor)
        #         self.vor_mids = np.asarray(vor_edge_mids) # edge midpoints

        self.n_vverts = len(self.vor_verts)
        self.n_vedges = len(self.vor_edges)
        self.n_vcells = len(self.vor_cells)

        # Calculate the centroid of the whole shape:
        self.centroid = np.mean(self.tri_verts, axis=0)

        # Define indices arrays for vor mesh:
        self.vor_vert_i = np.linspace(0, self.n_vverts - 1, self.n_vverts, dtype=np.int)
        self.vor_edge_i = np.linspace(0, self.n_vedges - 1, self.n_vedges, dtype=np.int)
        self.vor_cell_i = np.linspace(0, self.n_vcells - 1, self.n_vcells, dtype=np.int)

    def update_metrics(self):
        """
        Updates Voronoi cell area, and tri and vor edge lengths.
        This is required for case of deforming mesh, when
        connectivity of the meshes are not changing.

        """

        vor_cents = []  # centroid of voronoi polygons
        vor_sa = []  # surface area of voronoi polygons
        #         vor_edge_mids = [] # mids of Voronoi cell edges (same as tri edge mids)
        tri_edge_mids = []  # mids of tri cell edges
        vor_edge_len = []  # length of vor_edge
        tri_edge_len = []  # length of tri_edge

        tri_tang = []  # tangent vectors to tri_edges
        vor_tang = []  # tangent vectors to vor_edges

        for pts in self.vcell_verts:
            # Calculate centroid and area of the voronoi polygon:
            cx, cy = self.poly_centroid(pts)
            vor_cents.append([cx, cy])

            sa = self.area(pts)
            vor_sa.append(sa)

        self.vor_cents = np.asarray(vor_cents)
        self.vor_sa = np.asarray(vor_sa)

        for tei, (ti, tj) in enumerate(self.tri_edges):
            # find points representing vor and tri edges:

            tpi = self.tri_verts[ti]
            tpj = self.tri_verts[tj]

            tan_t = tpj - tpi

            tri_len = np.linalg.norm(tan_t)
            #     tri_len = np.sqrt(tan_t[0]**2 + tan_t[1]**2)

            pptm = (tpi + tpj) / 2  # calculate midpoint of tri-edge:

            tri_edge_mids.append(pptm)

            assert (tri_len != 0.0), "Tri-edge length equal to zero!"

            tri_edge_len.append(tri_len * 1)
            tri_tang.append(1 * tan_t / tri_len)

            # get the correct Voronoi edges
            # (they're at the same edge index as tri edges):

            vi, vj = self.vor_edges[tei]

            # get x-y points:
            vpi = self.vor_verts[vi]
            vpj = self.vor_verts[vj]

            tan_v = vpj - vpi

            vor_len = np.linalg.norm(tan_v)
            #     vor_len = np.sqrt(tan_v[0]**2 + tan_v[1]**2)

            #             ppvm = (vpi + vpj)/2 # midpoint of voronoi edge
            #             vor_edge_mids.append(ppvm)

            assert (vor_len != 0.0), "Vor-edge length equal to zero!"

            vor_edge_len.append(vor_len * 1)
            vor_tang.append(tan_v / vor_len)

        self.tri_tang = np.asarray(tri_tang)
        self.vor_tang = np.asarray(vor_tang)

        self.vor_edge_len = np.asarray(vor_edge_len)  # length of vor_edge
        self.tri_edge_len = np.asarray(tri_edge_len)  # length of tri_edge

        self.tri_mids = np.asarray(tri_edge_mids)  # edge midpoints
        #         self.vor_mids = np.asarray(vor_edge_mids) # edge midpoints

        # Calculate the centroid of the whole shape:
        self.centroid = np.mean(self.tri_verts, axis=0)

    def create_core_operators(self):
        """

        Creates the exterior derivative operators 'delta_0' and 'delta_1'.
        Note that the transpose of these matrices are equal to the
        boundary operators, where bount_1 = (delta_0).T and bound_2 = (delta_1).T.

        """

        # exterior derivative operator for tri mesh: operates on verts to return edges:
        delta_tri_0 = np.zeros((self.n_tedges, self.n_tverts))

        for ei, (vi, vj) in enumerate(self.tri_edges):
            delta_tri_0[ei, vj] = 1.0
            delta_tri_0[ei, vi] = -1.0

        self.delta_tri_0 = np.asarray(delta_tri_0)

        # get and store inverse:
        self.delta_tri_0_inv = np.linalg.pinv(self.delta_tri_0)

        # Mixing matrix inverse -- maps quantity from tri_edge mids to tri_verts using barycentric coordinates:
        # Note that forward mapping from tri_verts to tri_edge mids is given by M = np.abs(delta_tri_0)*(1/2),
        # therefore:
        self.MM_tri_inv = np.linalg.pinv(np.abs(self.delta_tri_0)*(1/2))

    def create_aux_operators(self):
        """
        Creates auxiliary operators required for curl, vector laplacians, etc. Note these are
        needed for the main mesh, but not for the 'mu-mesh' that is required for tensor work...

        """

        # exterior derivative operator for tri mesh operating on edges to return faces:
        delta_tri_1 = np.zeros((self.n_tcell, self.n_tedges))

        tri_edges = self.tri_edges.tolist()

        for ic, (vi, vj, vk) in enumerate(self.tri_cells):

            if [vi, vj] not in tri_edges:
                # get the index of the opposite sign edge:
                ea = tri_edges.index([vj, vi])
                delta_tri_1[ic, ea] = -1
            else:
                # get the forward sign edge:
                ea = tri_edges.index([vi, vj])
                delta_tri_1[ic, ea] = 1

            if [vj, vk] not in tri_edges:
                # get the index of the opposite sign edge:
                eb = tri_edges.index([vk, vj])
                delta_tri_1[ic, eb] = -1
            else:
                # get the forward sign edge:
                eb = tri_edges.index([vj, vk])
                delta_tri_1[ic, eb] = 1

            if [vk, vi] not in tri_edges:
                # get the index of the opposite sign edge:
                ec = tri_edges.index([vi, vk])
                delta_tri_1[ic, ec] = -1
            else:
                # get the forward sign edge:
                ec = tri_edges.index([vk, vi])
                delta_tri_1[ic, ec] = 1

        self.delta_tri_1 = np.asarray(delta_tri_1)

        # get and store inverse:
        self.delta_tri_1_inv = np.linalg.pinv(self.delta_tri_1)

        # exterior derivative operators for vor mesh: operates on verts to return edges
        delta_vor_0 = np.zeros((self.n_vedges, self.n_vverts))

        for ei, (vi, vj) in enumerate(self.vor_edges):
            delta_vor_0[ei, vj] = 1.0
            delta_vor_0[ei, vi] = -1.0

        self.delta_vor_0 = np.asarray(delta_vor_0)

        # get and store inverse:
        self.delta_vor_0_inv = np.linalg.pinv(self.delta_vor_0)

        ## Note the following creates the delta_vor_1 exterior derivative, however, delta_vor_1 is
        # not required as: delta_vor_1 = delta_tri_1.T. Keeping only for testing purposes.
        ## exterior derivative operator for vor mesh operating on edges to return faces:
        # delta_vor_1 = np.zeros((self.n_vcells, self.n_vedges))
        #
        # vor_edges = self.vor_edges.tolist()
        #
        # for ic, cell_verts in enumerate(self.vor_cells):
        #
        #     cell_verts_roll = np.roll(cell_verts, 1)
        #
        #     for (vi, vj) in zip(cell_verts, cell_verts_roll):
        #
        #         if [vi, vj] in vor_edges:
        #             # get the index of the opposite sign edge:
        #             ea = vor_edges.index([vi, vj])
        #             delta_vor_1[ic, ea] = 1  # Check which sign these should be depending on desired relations!
        #
        #         elif [vj, vi] in vor_edges:
        #             # get the forward sign edge:
        #             ea = vor_edges.index([vj, vi])
        #             delta_vor_1[ic, ea] = -1
        #
        # self.delta_vor_1 = np.asarray(delta_vor_1)

        self.MM_vor_inv = np.linalg.pinv(np.abs(self.delta_vor_0)*(1/2))

        # # Create mapping from tri verts to tri centers (Uses Barycentric coordinates to interpolate
        # # from verts to circumcentre):
        #
        # M_verts_to_cents = np.zeros((self.n_tcell, self.n_tverts))
        #
        # for ii, edge_inds in enumerate(self.tface_to_tedges):
        #     a, b, c = self.tri_edge_len[edge_inds]
        #
        #     b1o = (a ** 2) * (-a ** 2 + b ** 2 + c ** 2)
        #     b2o = (b ** 2) * (a ** 2 - b ** 2 + c ** 2)
        #     b3o = (c ** 2) * (a ** 2 + b ** 2 - c ** 2)
        #
        #     sumb = b1o + b2o + b3o
        #
        #     b1 = b1o / sumb
        #     b2 = b2o / sumb
        #     b3 = b3o / sumb
        #
        #     # get verts of triangle:
        #     vi, vj, vk = self.tri_cells[ii]
        #
        #     M_verts_to_cents[ii, vi] = b1
        #     M_verts_to_cents[ii, vj] = b2
        #     M_verts_to_cents[ii, vk] = b3
        #
        # self.M_verts_to_cents = np.asarray(M_verts_to_cents)
        # self.M_verts_to_cents_inv = np.linalg.pinv(self.M_verts_to_cents)


    #----Mathematical operator functions-----------

    def gradient_xy(self, S, gtype = 'tri'):
        """
        Gradient of scalar quantity 'S' with respect to the
        tangent vectors of tri_mesh (gtype = 'tri') or vor_mesh
        (gtype = 'vor').

        Note that this discrete gradient is a directional derivative with
        respect to the tangents of the mesh, and is not a true gradient in the x- and y-
        coordinate system.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if gradient is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        gradSx, gradSy  -- the x and y components of the directional derivative of S

        """

        if gtype == 'tri':

            assert(len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            gradSx, gradSy = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, S)*self.tri_tang.T

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            assert(len(S) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            gradSx, gradSy = (1/self.vor_edge_len)*np.dot(self.delta_vor_0, S)*self.vor_tang.T

        else:
            gradSx = None
            gradSy = None


        return gradSx, gradSy

    def gradient(self, S, gtype = 'tri'):
        """
        Gradient of scalar quantity 'S' with respect to the
        tangent vectors of tri_mesh (gtype = 'tri') or vor_mesh
        (gtype = 'vor').

        Note that this discrete gradient is a directional derivative with
        respect to the tangents of the mesh, and returns the component of the
        gradient with respect to the tangent vectors of the mesh on which it was
        computed.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if gradient is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        gradS  -- the tangential component of the directional derivative of S along tangents of mesh

        """

        if gtype == 'tri':

            assert(len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            gradS = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, S)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            assert(len(S) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            gradS = (1/self.vor_edge_len)*np.dot(self.delta_vor_0, S)

        else:
            gradS = None

        return gradS

    def div_xy(self, Fx, Fy, gtype = 'tri'):
        """
        Divergence of a vector field Fx, and Fy with respect to the tri (gtype = 'tri') or vor (gtype = 'vor')
        mesh cells.

        If gtype == 'tri', the divergence is with respect to a vor_mesh "control cell", and the final
        divergence is defined at vor_cents (tri vertices).

        If gtype == 'vor', the divergence is with respect to a tri_mesh "control cell", and the final
        divergence is defined at tri_cents (vor vertices).

        This computation assumes that self.delta_vor_1 = -self.delta_tri_0.T, which has been confirmed,
        and can be confirmed by uncommenting code calculating self.delta_vor_1 in create_aux_operators method, above.

        Parameters
        -----------
        Fx, Fy   -- components of vector field defined on mesh edges, depending on gtype
        gtype -- specifies if divergence is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        divF  -- divergence of the vector field Fx, Fy

        """

        if gtype == 'tri':

            # get component of vector field parallel to tri_mesh tangents (and therefore perpendicular to vor_edges):
            FF = Fx*self.tri_tang[:,0] + Fy*self.tri_tang[:,1]

            divF = (1/self.vor_sa)*np.dot(-self.delta_tri_0.T, self.vor_edge_len*FF)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor div"

            # get component of vector field parallel to vor_mesh tangents (and therefore perpendicular to tri_edges):
            FF = Fx*self.vor_tang[:,0] + Fy*self.vor_tang[:,1]

            divF = (1/self.tri_sa)*np.dot(self.delta_tri_1, self.tri_edge_len*FF)

        else:
            divF = None


        return divF

    def div(self, Ft, gtype = 'tri'):
        """
        Divergence of a vector field tangential component Ft with respect to the tri (gtype = 'tri')
        or vor (gtype = 'vor') mesh cells.

        If gtype == 'tri', the divergence is with respect to a vor_mesh "control cell", and the final
        divergence is defined at vor_cents (tri vertices). The input Ft should represent the tangential component
        of a vector field with respect ot the tri_tangents.

        If gtype == 'vor', the divergence is with respect to a tri_mesh "control cell", and the final
        divergence is defined at tri_cents (vor vertices). The input Ft should represent the tangential component
        of a vector field with respect ot the vor_tangents.

        This computation assumes that self.delta_vor_1 = -self.delta_tri_0.T, which has been confirmed,
        and can be confirmed by uncommenting code calculating self.delta_vor_1 in create_aux_operators method, above.

        Parameters
        -----------
        Ft   -- tangential component of vector field with respect to mesh edges, depending on gtype
        gtype -- specifies if divergence is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        divF  -- divergence of the vector field Ft

        """

        if gtype == 'tri':

            divF = (1/self.vor_sa)*np.dot(-self.delta_tri_0.T, self.vor_edge_len*Ft)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor div!"

            divF = (1/self.tri_sa)*np.dot(self.delta_tri_1, self.tri_edge_len*Ft)

        else:
            divF = None


        return divF

    def lap(self, S, gtype = 'tri'):
        """
        Computes a scalar forwards Laplacian on a scalar variable S as the divergence of the gradient of the
        scalar.

        If gtype = 'tri', the gradient is taken with respect to the tri_mesh edges with vor_mesh control volumes,
        and if gtype = 'vor', the gradient is taken with respect to the vor_mesh edges with tri_mesh control volumes.

        Note that due to the structure of the grids, the gtype tri lap is closed boundary (zero flux) while the
        gtype vor lap is open boundary.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if laplacian is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        lapS  -- the Laplacian of S with 'natural' boundary conditions.

        """
        if gtype == 'tri':

            # calculate gradient of S:
            gS = self.gradient(S, gtype='tri')

            # calculate the divergence of the gradient, which is the laplacian:
            lapS = self.div(gS, gtype = 'tri')

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate gradient of S:
            gS = self.gradient(S, gtype='vor')

            # calculate the divergence of the gradient, which is the Laplacian:
            lapS = self.div(gS, gtype = 'vor')

        else:
            lapS = None

        return lapS

    def lap_inv(self, S, gtype = 'tri'):
        pass

    def curl(self, Fz, gtype = 'tri'):
        pass

    def curl_xy(self, Fx, Fy, gtype = 'tri'):
        pass

    def verts_to_mids(self, Sv, gtype = 'tri'):
        """
        Maps property S from vertices of mesh to midpoints of edges using barycentric interpolation.

        Parameters
        ------------
        Sv -- property defined at vertices (verts depend on gtype)
        gtype  -- if transformation is from triverts to trimids or vorverts to vormids

        Returns
        --------
        Sm  -- property at edge mids

        """

        if gtype == 'tri':
            assert(len(Sv) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"
            MM = np.abs(self.delta_tri_0)*(1/2)

            Sm = np.dot(MM, Sv)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"
            assert(len(Sv) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            MM = np.abs(self.delta_vor_0)*(1/2)

            Sm = np.dot(MM, Sv)

        else:
            Sm = None

        return Sm

    def mids_to_verts(self, Sm, gtype = 'tri'):
        """
        Maps property Sm from edge mids of mesh to vertices using barycentric interpolation.

        Parameters
        ------------
        Sm -- property at edge mids
        gtype  -- if transformation is from triverts to trimids or vorverts to vormids

        Returns
        --------
        Sv  -- property defined at vertices (verts depend on gtype)

        """

        if gtype == 'tri':
            assert(len(Sm) == self.n_tedges), "Length of array passed to gradient is not edges length"
            MM_inv = self.MM_tri_inv

            Sv = np.dot(MM_inv, Sm)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"
            assert(len(Sm) == self.n_vedges), "Length of array passed to gradient is not edges length"

            MM_inv = self.MM_vor_inv

            Sv = np.dot(MM_inv, Sm)

        else:

            Sv = None

        return Sv

    def verts_to_cent(self, Sv):

        pass

    def cent_to_verts(self, Sc):

        pass

    def vector_laplacian(self, Fx, Fy, gtype = 'tri'):
        pass

    def vector_laplacian_inv(self, Fx, Fy, gtype = 'tri'):
        pass

    def helmholtz_hodge(self, Fx, Fy, gtype = 'tri'):

        pass

    def calc_tri(self):
        self.tri_circcents = []
        self.tri_rcircs = [] # circumradius of triangle
        self.tri_rin = [] # inradius of triangle
        self.tri_sa = []  # surface area of triangle faces

        for i, vert_inds in enumerate(self.tri_cells):

            abc = self.tri_verts[vert_inds]
            vx, vy, r_circ, r_in = self.circumc(abc[0], abc[1], abc[2])

            sa = self.area(abc)  # surface area of triangle

            self.tri_circcents.append([vx, vy])
            self.tri_rcircs.append(r_circ)  # circumradius of triangle
            self.tri_rin.append(r_in)  # inradius of triangle
            self.tri_sa.append(sa)  # surface area of triangle faces

        self.tri_circcents = np.asarray(self.tri_circcents)
        self.tri_rcircs = np.asarray(self.tri_rcircs) # circumradius of triangle
        self.tri_rin = np.asarray(self.tri_rin) # inradius of triangle
        self.tri_sa = np.asarray(self.tri_sa)  # surface area of triangle faces

    def mesh_quality_calc(self):

        uu = (1/2)*self.tri_rcircs*(self.tri_rcircs - 2*self.tri_rin)

        # energy at edge cells is double that of those in the interior:
        uu[self.bflags_tcells] = self.tri_rcircs[self.bflags_tcells]*(self.tri_rcircs[self.bflags_tcells]
                                                                      - 2*self.tri_rin[self.bflags_tcells])

        return uu

    def search_point_cloud(self, pts, pt_cloud):

        if len(pts) == 0:

            matched_pts = []
            unmatched_pts = []

        elif len(pt_cloud) == 0:

            matched_pts = []
            unmatched_pts = []

        else:

            search_M = 99 * np.ones((len(pts), len(pt_cloud)))

            pts_inds = np.asarray([i for i, cc in enumerate(pts)])

            for i, pt in enumerate(pts):

                for j, ptc in enumerate(pt_cloud):

                    dist = np.sqrt((pt[0] - ptc[0]) ** 2 + (pt[1] - ptc[1]) ** 2)

                    if dist < 1.0e-15:

                        search_M[i, j] = 1.0

                    else:

                        search_M[i, j] = 0.0

            matched_inds, matched_inds_cloud = (search_M == 1).nonzero()

            unmatched_inds = np.setdiff1d(pts_inds, matched_inds)

            matched_pts = pts[matched_inds]
            unmatched_pts = pts[unmatched_inds]

        return matched_pts, unmatched_pts

    def refine_mesh(self, max_steps=25, convergence=7.5):

        logs.log_info("Initializing Voronoi mesh optimization...")

        opti_steps = np.arange(max_steps)

        ui = self.mesh_quality_calc()

        UU = np.sum(ui)/self.cell_radius**2

        for i in opti_steps:

            if UU > convergence:

                # Continuously reassign tri_verts to vor_centres, without affecting the boundary
                self.tri_verts[self.biocell_i] = self.vor_cents[self.biocell_i]

                self.create_tri_mesh()
                self.create_mappings()
                self.process_vormesh()

                ui = self.mesh_quality_calc()

                UU = np.sum(ui)/self.cell_radius**2

                conv_mess = "Step {}: mesh quality {}".format(i, UU)
                #                 logs.log_info(conv_mess)
                logs.log_info(conv_mess)

            else:

                # Finish up:
                self.init_mesh() # build entire mesh

                self.mesh_qual = UU
                #                 logs.log_info("Convergence condition met for mesh optimization.")
                print("Convergence condition met for mesh optimization.")
                final_mess = "Final mesh quality {}".format(UU)
                #                 logs.log_info(final_mess)
                logs.log_info(final_mess)
                break

    def clip_to_curve(self, imagemask):

        """
        Uses an image mask clipping curve and Sutherland Hodgmann algorithm to clip the
        cells of a voronoi mask at the curve boundary. It then recalculates the mesh with the
        voronoi cell centers used as trimesh verts. A good algorithm, but only works for
        concave clipping curves.

        :param imagemask:
        :return:
        """

        clip_vor_verts = []
        clip_vor_cents = []

        for ii, (poly_ind, cell_poly, vor_cent) in enumerate(zip(self.vor_cells,
                                                                 self.vcell_verts,
                                                                 self.vor_cents)):

            if len(poly_ind) >= 3:
                cell_polya = cell_poly.tolist()
                point_check = np.zeros(len(cell_poly))

                for i, pnt in enumerate(cell_poly):

                    point_val = imagemask.clipping_function(pnt[0], pnt[1])

                    if point_val != 0.0:
                        point_check[i] = 1.0

                if point_check.sum() == len(cell_poly):  # if all points are all inside the clipping zone

                    clip_vor_verts.append(np.array(cell_polya))
                    cx, cy = self.poly_centroid(cell_polya)
                    clip_vor_cents.append([cx, cy])

                elif point_check.sum() > 0.0 and point_check.sum() < len(
                        cell_poly):  # the region's points are in the clipping func range

                    clip_poly = clip_counterclockwise(
                        cell_poly, imagemask.clipcurve)

                    if len(clip_poly):

                        # For inside voronoi cell: augment old_poly_pts with keeper_pts
                        old_poly_pts, new_poly_pts = self.search_point_cloud(clip_poly, cell_poly)

                        throw_away_pts, keeper_pts = self.search_point_cloud(new_poly_pts, imagemask.clipcurve)

                        # For outside voronoi cell: augment outside_cell_b with keeper_pts
                        outside_cell_a, outside_cell_b = self.search_point_cloud(cell_poly, clip_poly)

                        if len(keeper_pts) == 0:

                            out_stack = []
                            in_stack = old_poly_pts * 1

                        else:

                            out_stack = np.vstack((outside_cell_b, keeper_pts))
                            in_stack = np.vstack((old_poly_pts, keeper_pts))

                        if len(in_stack) >= 3:
                            in_vor_cell = orient_counterclockwise(in_stack)
                            clip_vor_verts.append(in_vor_cell)

                            cx, cy = self.poly_centroid(in_vor_cell)
                            clip_vor_cents.append([cx, cy])


        # clip_vor_verts = np.asarray(clip_vor_verts)
        clip_vor_cents = np.asarray(clip_vor_cents)

        self.tri_verts = clip_vor_cents
        self.image_mask = imagemask

        self.init_mesh()


    #-----Utility functions--------------------------

    def area(self, p):
        """
        Calculates the unsigned area of an arbitrarily shaped polygon defined by a set of
        counter-clockwise oriented points in 2D.

        Parameters
        ----------
        p               xy list of polygon points

        Returns
        -------
        area            area of a polygon in square meters

        Notes
        -------
        The algorithm is an application of Green's theorem for the functions -y and
        x, exactly in the way a planimeter works.
        """

        foo = np.asarray(p)

        # move points along by one:
        foo_p = np.roll(foo, -1, axis=0)

        ai = foo[:, 0] * foo_p[:, 1] - foo_p[:, 0] * foo[:, 1]

        # aa = np.abs((1 / 2) * np.sum(ai))

        aa = (1 / 2) * np.sum(ai)

        return aa

    def poly_centroid(self, p):
        """
        Calculates the centroid (geometric centre of mass) of a polygon.

        Parameters
        ----------
        p       array of [x,y] points defining polygon vertices

        Returns
        --------
        cx, cy  polygon centroid coordinates

        reference: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
        """

        foo = np.asarray(p)

        # move points along by one:
        foo_p = np.roll(foo, -1, axis=0)

        ai = foo[:, 0] * foo_p[:, 1] - foo_p[:, 0] * foo[:, 1]

        aa = (1 / 2) * np.sum(ai)  # signed area

        cx = (1 / (6 * aa)) * np.sum((foo[:, 0] + foo_p[:, 0]) * (foo[:, 0] * foo_p[:, 1] - foo_p[:, 0] * foo[:, 1]))
        cy = (1 / (6 * aa)) * np.sum((foo[:, 1] + foo_p[:, 1]) * (foo[:, 0] * foo_p[:, 1] - foo_p[:, 0] * foo[:, 1]))

        return cx, cy


    def circumc(self, A, B, C):
        """
        Calculates the circumcenter and circumradius of a triangle with
        vertices A = [Ax, Ay], B = [Bx, By], and C = [Cx, Cy]

        returns ox, oy, rc, the x and y coordinates of the circumcentre and
        the circumradius, respectively.

        """

        # Point coords:
        Ax = A[0]
        Ay = A[1]
        Bx = B[0]
        By = B[1]
        Cx = C[0]
        Cy = C[1]

        # Calculate circumcentre:
        # (from https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_2)

        A2 = Ax ** 2 + Ay ** 2
        B2 = Bx ** 2 + By ** 2
        C2 = Cx ** 2 + Cy ** 2

        denom = 2 * (Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))

        ox = (A2 * (By - Cy) + B2 * (Cy - Ay) + C2 * (Ay - By)) / denom
        oy = (A2 * (Cx - Bx) + B2 * (Ax - Cx) + C2 * (Bx - Ax)) / denom

        # Calculate circumradius:
        # (from https://www.mathalino.com/reviewer/
        # derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle)
        a = np.sqrt((Ax - Bx) ** 2 + (Ay - By) ** 2)
        b = np.sqrt((Bx - Cx) ** 2 + (By - Cy) ** 2)
        c = np.sqrt((Cx - Ax) ** 2 + (Cy - Ay) ** 2)

        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))

        # circumcircle:
        if area > 0.0:
            rc = a * b * c / (4.0 * area)

        else:
            rc = 0.0

        # inradius:
        ri = (2*area)/ (a + b + c)

        return ox, oy, rc, ri


    #WasteLands-------------------------------------
    # # Interpolate mesh energy from tri cents to verts
    # uv = np.dot(self.M_verts_to_cents_inv, ui)
    # # Calculate the gradient (x and y components on tri edges):
    # guxe, guye = self.gradient_xy(uv)
    # # Interpolate from edges to verts:
    # gUx = self.mids_to_verts(guxe)
    # gUy = self.mids_to_verts(guye)
    #
    # # move the verts to the negative gradient:
    # self.tri_verts[self.biocell_i,0] += -gUx[self.biocell_i]*step_size
    # self.tri_verts[self.biocell_i,1] += -gUy[self.biocell_i]*step_size
    #
    # # recalculate key features of tri_mesh with new triverts positions
    # # self.calc_tri()




