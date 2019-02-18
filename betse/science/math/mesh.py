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

from betse.science.math.geometry.polygon.geopolyconvex import clip_counterclockwise
from betse.science.math.geometry.polygon.geopoly import orient_counterclockwise

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
                 image_mask=None):

        self.single_cell_noise = single_cell_noise
        self.single_cell_sides = single_cell_sides
        self.alpha_shape = alpha_shape
        self.cell_radius = cell_radius
        self.use_alpha_shape = use_alpha_shape
        self.image_mask = image_mask

        if seed_points is None:
            self.make_single_cell_points()
        else:
            self.tri_verts = seed_points


        self.init_mesh()

    def init_mesh(self):

        self.create_tri_mesh()
        self.create_mappings()
        self.process_vormesh()
        self.create_d_operators()

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

        # Calculate the centroids of the triangulation (which should be identical to mem_verts,
        # as well as circumcircle radii:

        self.tri_circcents = []
        self.vor_verts_bound = []
        self.tri_rcircs = []
        self.tricell_i = []  # index of simplexes
        self.tri_sa = []  # surface area of triangle faces
        self.tcell_verts = []  # x,y components of tri_mesh cells
        simplices2 = []

        for i, vert_inds in enumerate(trimesh.simplices):

            abc = trimesh.points[vert_inds]
            vx, vy, r_circ = self.circumc(abc[0], abc[1], abc[2])

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
                        self.tricell_i.append(i)
                        simplices2.append(vert_inds)
                        self.tri_sa.append(sa)
                        self.tcell_verts.append(abc)

                self.tri_cells = np.asarray(simplices2)  # reassign point inds to retained simplices

            else:
                self.tri_cells = trimesh.simplices
                self.tri_circcents.append([vx, vy])
                self.tri_rcircs.append(r_circ)
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

        self.tface_to_tedges = np.asarray(face_to_edges) # tri_face index to tri_edges indices mapping

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

            assert (len(shared_ij) == 2), "Shared vor_cell edge inds must be length 2"

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

    def create_d_operators(self):
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

        # Mixing matrix inverse -- maps quantity from tri_edge mids to tri_verts using barycentric coordinates:
        # Note that forward mapping from tri_verts to tri_edge mids is given by M = np.abs(delta_tri_0)*(1/2),
        # therefore:
        self.MM_tri_inv = np.linalg.pinv(np.abs(self.delta_tri_0)*(1/2))
        self.MM_vor_inv = np.linalg.pinv(np.abs(self.delta_vor_0)*(1/2))


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

        aa = np.abs((1 / 2) * np.sum(ai))

        # aa = 0.5 * abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1, y1)) in zip(p, p[1:] + [p[0]])))

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

        return ox, oy, rc

    def mesh_energy_i(self, A, B, C, ro = 5.0e-5):
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

        # Calculate circumradius:
        # (from https://www.mathalino.com/reviewer/
        # derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle)
        a = np.sqrt((Ax - Bx) ** 2 + (Ay - By) ** 2)
        b = np.sqrt((Bx - Cx) ** 2 + (By - Cy) ** 2)
        c = np.sqrt((Cx - Ax) ** 2 + (Cy - Ay) ** 2)

        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))

        rt = (1/ro)*2*area/(a + b + c)

        Rt = (1/ro)*(a*b*c)/(4*area)

        if area < 0.0:

            uu = 1.0e10

        else:

            uu = (1/2)*Rt*(Rt - 2*rt)


        return uu

    def mesh_quality_calc(self):

        edge_lengths = self.vor_edge_len
        mean_edge = edge_lengths.mean()

        edge_qual = (edge_lengths / mean_edge) * 100
        mesh_qual_metric = edge_qual.min()

        Ui = []

        for ci, (vi, vj, vk) in enumerate(self.tri_cells):
            # get points
            pi, pj, pk = self.tri_verts[[vi, vj, vk]]

            ui = self.mesh_energy_i(pi, pj, pk)

            Ui.append(ui)
        #
        UU = np.sum(Ui)

        return mesh_qual_metric, UU

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

    def refine_mesh(self, max_steps=50, convergence=15.0):

        #         logs.log_info("Initializing Voronoi mesh optimization...")

        # optimization vector:
        #         opti_steps = np.arange(p.maximum_voronoi_steps)

        opti_steps = np.arange(max_steps)

        mesh_qual_metric, UU = self.mesh_quality_calc()

        # fixed_bound_points = self.tri_mids[self.bflags_tedges] # fix a bounding curve of points

        # self.tri_verts = np.vstack((fixed_bound_points, self.vor_cents)) # add in the fixed-bound points once!

        # # change the points at the boundary by folding in tri_mids, just once:
        # self.tri_verts = np.vstack((self.tri_mids[self.bflags_tedges], self.tri_verts)) # add in the fixed-bound points
        # self.create_tri_mesh()
        # self.create_mappings()
        # self.process_vormesh()

        for i in opti_steps:

            if mesh_qual_metric < convergence:

                # Continuously reassign tri_verts to vor_centres, without affecting the boundary
                self.tri_verts[self.biocell_i] = self.vor_cents[self.biocell_i]

                self.create_tri_mesh()
                self.create_mappings()
                self.process_vormesh()

                mesh_qual_metric, UU = self.mesh_quality_calc()

                conv_mess = "Step {}: mesh quality {}, {}".format(i, mesh_qual_metric, UU)
                #                 logs.log_info(conv_mess)
                print(conv_mess)

            else:

                # Finish up:
                self.create_d_operators() # calculate the matrices

                self.mesh_qual = mesh_qual_metric
                #                 logs.log_info("Convergence condition met for mesh optimization.")
                print("Convergence condition met for mesh optimization.")
                final_mess = "Final mesh quality {}".format(mesh_qual_metric)
                #                 logs.log_info(final_mess)
                print(final_mess)
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


### WASTELANDS ###

    #         # exterior derivative operator for vor mesh operating on edges to return faces:
    #         delta_vor_1 = np.zeros((self.n_vcells, self.n_vedges))

    #         vor_edges = self.vor_edges.tolist()

    #         for ic, cell_verts in enumerate(self.vor_cells):

    #             cell_verts_roll = np.roll(cell_verts, 1)

    #             for (vi, vj) in zip(cell_verts, cell_verts_roll):

    #                 if [vi, vj] in vor_edges:
    #                     # get the index of the opposite sign edge:
    #                     ea = vor_edges.index([vi, vj])
    #                     delta_vor_1[ic, ea] = -1

    #                 elif [vj, vi] in vor_edges:
    #                     # get the forward sign edge:
    #                     ea = vor_edges.index([vj, vi])
    #                     delta_vor_1[ic, ea] = 1

    #         self.delta_vor_1 = np.asarray(delta_vor_1)