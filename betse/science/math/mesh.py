#!/usr/bin/env python3
# ....................{ LICENSE                           }....................
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................

import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import colors
from matplotlib import colorbar
from matplotlib import rcParams
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Circle
from matplotlib import path

from betse.util.math.geometry.polygon.geopolyconvex import clip_counterclockwise
from betse.util.math.geometry.polygon.geopoly import orient_counterclockwise, is_convex, is_cyclic_quad
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
                 cell_radius = None, # average half distance between seed points
                 mesh_type = 'tri', # mesh_type can be 'tri' or 'quad'
                 seed_points=None, # points to use in Delaunay triangulation
                 use_centroids = True, # Use the centroids of polygons as Voronoi verts instead of circumcenters
                 use_alpha_shape=True, # use alpha-shape exclusion of triangles from Delaunay
                 alpha_shape=0.4, # Threshhold for alpha-shape analysis
                 single_cell_sides=6,  # If seed points is None, a single cell will be created with this many sides
                 single_cell_noise=0.5,  # Irregularity to add to single cell construction
                 image_mask=None, # Use an image mask (from BETSE) to control point location
                 make_all_operators = True, # Make all operators (True), or only key ones (False)
                 allow_merging = True, # Allow tri-cells to be merged to quads if circumcenters are close?
                 merge_thresh = 0.2, # Distance threshhold (%of total radius) for merging close circumcenters
                 close_thresh = 0.25, # threshhold for removal of close tri vert neighour points
                 ):


        self.single_cell_noise = single_cell_noise
        self.single_cell_sides = single_cell_sides
        self.alpha_shape = alpha_shape
        self.use_alpha_shape = use_alpha_shape
        self.image_mask = image_mask
        self.make_all_operators = make_all_operators
        self.mesh_type = mesh_type
        self.merge_thresh = merge_thresh
        self.allow_merging = allow_merging
        self.use_centroids = use_centroids
        self.close_thresh = close_thresh

        self.removed_bad_verts = False # flag to bad tri-vert removal (only do case once!)

        if cell_radius is None and seed_points is not None:
            # find the average distance between points in the cluster:
            point_tree = cKDTree(seed_points)
            nd, _ = point_tree.query(seed_points, k=2)
            self.cell_radius = np.mean(nd[:,1])/2

        elif cell_radius is None and seed_points is None:
            self.cell_radius = 5.0e-6

        else:
            self.cell_radius = cell_radius

        if seed_points is None:
            self.make_single_cell_points()
        else:
            self.tri_verts = seed_points

    def init_mesh(self):

        self.create_tri_mesh()
        self.create_mappings()
        self.process_voredges()
        self.process_vorcells()
        self.create_core_operators()

        if self.make_all_operators:
            self.create_aux_operators()

    def pre_mesh(self):
        """
        Create only the tri and vor meshes without creating operators.

        """
        self.create_tri_mesh()
        self.create_mappings()
        self.process_voredges()
        self.process_vorcells()

    def init_and_refine(self, smoothing = None, refinement = True,
                        max_steps=25, convergence=7.5, fix_bounds=True):

        if smoothing is not None:
            self.init_mesh()  # init the whole mesh
            self.laplacian_smoothing(stepsize=smoothing)  # run the smoothing of the tri_verts
            self.pre_mesh()

        else:
            self.pre_mesh()

        if refinement:

            self.refine_mesh(max_steps=max_steps, convergence=convergence, fix_bounds=fix_bounds)
            self.pre_mesh()

        self.create_core_operators()

        if self.make_all_operators:
            self.create_aux_operators()

    def clip_and_refine(self, imagemask, smoothing = None, refinement = True,
                        max_steps=25, convergence=7.5, fix_bounds=True):

        self.pre_mesh()
        self.clip_to_curve(imagemask)

        if smoothing is not None:
            self.init_mesh()  # init the whole mesh to remake operators
            self.laplacian_smoothing(stepsize=smoothing)  # run the smoothing of the tri_verts
            self.pre_mesh()

        if refinement:

            self.refine_mesh(max_steps=max_steps, convergence=convergence, fix_bounds=fix_bounds)
            self.pre_mesh()

        self.create_core_operators()

        if self.make_all_operators:
            self.create_aux_operators()

    def make_single_cell_points(self):
        """
        Creates seed points (tri_mesh verts) for a single closed voronoi cell and
        its boundary.

        """

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

        self.trimesh_core_calcs()

        self.sanity_check() # check for bad triverts (triverts belonging to no simplexes)

        # Next check for triverts with really close circumcenters, and, if necessary, merge to quads,
        # or if requested, try merging as much of the mesh to quads as possible
        # Find ccents that are close to one another and mark cells for quad-merge
        tri_tree = cKDTree(self.tri_ccents)
        di, ni = tri_tree.query(self.tri_ccents, k=2)
        mark_for_merge = []
        for si, ccenti in enumerate(self.tri_ccents):
            disti = di[si, 1] # distance between circumcenter of si and nearest neighbour
            indi = ni[si, 1] # index to nearest neighbour of si

            if self.mesh_type == 'quad': # if user wants a quad mesh anyway, append all
                # nearest-neighbour simplices:
                mark_for_merge.append([indi])

            else:
                # only mark if distance between ccents is less than 20% of the cell radius
                if disti < self.cell_radius * self.merge_thresh:
                    mark_for_merge.append([indi])

                else: # append an empty list to indicate there's no merging to be performed
                    mark_for_merge.append([])

        if self.allow_merging or self.mesh_type == 'quad':
            # merge tri elements as required:
            self.merge_tri_mesh(mark_for_merge)

        # process the edges and boundaries:
        self.process_primary_edges(mesh_version = self.mesh_type)

        # inds to full voronoi cells only (these are the official cells of cluster):
        self.biocell_i = np.delete(self.tri_vert_i, self.bflags_tverts)
        # Calculate the centroid of the whole shape:
        self.centroid = np.mean(self.tri_verts, axis=0)

    def trimesh_core_calcs(self):

        # Figure out if there are super close points in the tri_verts set:
        # Next check for triverts with really close circumcenters, and, if necessary, merge to quads,
        # or if requested, try merging as much of the mesh to quads as possible
        # Find ccents that are close to one another and mark cells for quad-merge
        tri_tree = cKDTree(self.tri_verts)
        di, ni = tri_tree.query(self.tri_verts, k=2)
        mark_for_merge = set()
        not_merged = set()
        for vi, tvertpts in enumerate(self.tri_verts):
            disti = di[vi, 1]  # distance between circumcenter of si and nearest neighbour
            indi = ni[vi, 1]  # index to nearest neighbour of si
            near_point_check = disti/self.cell_radius
            if near_point_check < self.close_thresh:
                if vi not in not_merged:
                    mark_for_merge.add(vi)
                not_merged.add(indi)

        mark_for_merge = np.asarray(list(mark_for_merge))
        self.tri_verts = np.delete(self.tri_verts, mark_for_merge, axis = 1)

        # calculate the Delaunday triangulation based on the cluster-masked seed points:
        trimesh = Delaunay(self.tri_verts)  # Delaunay trianulation of cell centres
        self.n_tverts = len(self.tri_verts)  # number of tri_verts
        self.tri_vert_i = np.linspace(0, self.n_tverts - 1,
                                      self.n_tverts, dtype=np.int)

        tri_ccents = []  # circumcentres of the triangles
        tri_cents = []  # centroids of the triangles
        tri_rcircs = []  # circumradius of triangle
        tri_sa = []  # surface area of triangle faces
        tcell_verts = []  # x,y components of tri_mesh cells
        tri_cells = []  # indices to tri_verts defining each triangle (simplex)

        for si, vert_inds in enumerate(trimesh.simplices):

            abc = trimesh.points[vert_inds]
            vx, vy, r_circ, r_in = self.circumc(abc[0], abc[1], abc[2])
            cx, cy = self.poly_centroid(abc)

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
                        tri_ccents.append([vx, vy])
                        tri_cents.append([cx, cy])
                        tri_rcircs.append(r_circ)
                        tri_cells.append(vert_inds)
                        tri_sa.append(sa)
                        tcell_verts.append(abc)

            else:  # append all information without screening for triangle suitability
                tri_ccents.append([vx, vy])
                tri_cents.append([cx, cy])
                tri_rcircs.append(r_circ)
                tri_sa.append(sa)
                tri_cells.append(vert_inds)
                tcell_verts.append(abc)

        self.tri_cells = np.asarray(tri_cells)  # reassign point inds to retained simplices

        self.n_tcell = len(tri_cells)  # number of simplexes in trimesh
        self.tri_cell_i = np.asarray([i for i in range(self.n_tcell)])  # indices vector of trimesh

        self.tri_ccents = np.asarray(tri_ccents)
        self.tri_cents = np.asarray(tri_cents)
        self.tri_rcircs = np.asarray(tri_rcircs)
        self.tri_sa = np.asarray(tri_sa)
        self.tcell_verts = np.asarray(tcell_verts)

        # Create an updated mapping of which triangle each vertex belongs to:
        self.create_tri_map()

    def merge_tri_mesh(self, merge_list):
        """
        Merge cells of the tri_mesh into a (potentially) mixed quad-tri mesh.

        """

        quad_cells = []
        quad_ccents = []
        quad_sa = [] # area of quadrilateral cell
        quad_rcircs = [] # circumradius of quad
        quadcell_i = []  # index of simplexes

        qcell_verts = []  # verts of tri or quad simplex
        quad_cents = [] # centroids of tri or quad simplex

        # for triangle index a, find triangles with nearly coincident circumcentres:
        for ai, blist in enumerate(merge_list):

            vertsa = self.tri_cells[ai] # verts of triangle ai
            cc1 = self.tri_ccents[ai] # ccent of triangle ai
            c1 = self.tri_cents[ai] # cent of triangle ai
            a1 = self.tri_sa[ai]  # surface area of triangle ai
            rc1 = self.tri_rcircs[ai]  # circumcenter of triangle ai

            if len(blist) != 0:

                bi = blist[0]

                vertsb = self.tri_cells[bi]

                cc2 = self.tri_ccents[bi]
                c2 = self.tri_cents[bi]
                a2 = self.tri_sa[bi]
                rc2 = self.tri_rcircs[bi]

                # If triangle ai has not yet been used in a merging:
                if ai in self.free_to_merge:

                    # triangle bi has also not yet been used in a merging:
                    if bi in self.free_to_merge:
                        # get verts of tri a and tri b:
                        quad_i = np.unique((vertsa, vertsb))

                        # also get the shared verts, which represent the shared edge:
                        shared_ij = np.intersect1d(vertsa, vertsb)

                        # if the resulting merger leads to 4 unique vertices and one edge:
                        if len(quad_i) == 4 and len(shared_ij) == 2:
                            # orient verts counterclockwise:
                            quad_pts = self.tri_verts[quad_i]
                            cent = quad_pts.mean(axis=0)  # calculate the centre point
                            angles = np.arctan2(quad_pts[:, 1] - cent[1],
                                                quad_pts[:, 0] - cent[0])  # calculate point angles
                            sorted_region = quad_i[np.argsort(angles)]  # sort indices counter-clockwise

                            sorted_pts = quad_pts[np.argsort(angles)]

                            # test to see if the merged poly is convex:
                            conv_quad = is_convex(sorted_pts) #
                            cycl_quad = is_cyclic_quad(sorted_pts[0], sorted_pts[1],
                                                                          sorted_pts[2], sorted_pts[3])

                            if conv_quad and cycl_quad:

                                Rq, areaq, ccxq, ccyq = self.quad_circumc(sorted_pts[0], sorted_pts[1],
                                                                          sorted_pts[2], sorted_pts[3])

                                quad_cells.append(sorted_region)
                                qcell_verts.append(sorted_pts)

                                quad_ccents.append([ccxq, ccyq])
                                quad_rcircs.append(Rq)

                                quad_sa.append(areaq)

                                quadcell_i.append(ai)

                                # Calculate centroid of the quad cell:
                                cx, cy = self.poly_centroid(sorted_pts)
                                quad_cents.append([cx, cy])


                            else: # if it's not convex and/or cyclic, then don't split it; append original triangles:

                                ptsa = self.tri_verts[vertsa] # get the x, y coords of points
                                ptsb = self.tri_verts[vertsb]

                                quad_cells.append(vertsa)
                                qcell_verts.append(ptsa)
                                quad_ccents.append(cc1)
                                quad_rcircs.append(rc1)
                                quad_sa.append(a1)
                                quadcell_i.append(ai)
                                quad_cents.append(c1)

                                quad_cells.append(vertsb)
                                qcell_verts.append(ptsb)
                                quad_ccents.append(cc2)
                                quad_rcircs.append(rc2)
                                quad_sa.append(a2)
                                quadcell_i.append(bi)
                                quad_cents.append(c2)

                            # remove triangles ai and bi from future merges:
                            self.free_to_merge.remove(ai)
                            self.free_to_merge.remove(bi)

                    else:  # if bi is not in free to merge, add in ai as a triangle:

                        ptsa = self.tri_verts[vertsa]  # get the x, y coords of points

                        quad_cells.append(vertsa)
                        qcell_verts.append(ptsa)
                        quad_ccents.append(cc1)
                        quad_rcircs.append(rc1)
                        quad_sa.append(a1)
                        quadcell_i.append(ai)
                        quad_cents.append(c1)

                        self.free_to_merge.remove(ai)

            else: # if there's no request to merge, append all original triangle features to the list

                ptsa = self.tri_verts[vertsa]  # get the x, y coords of points

                quad_cells.append(vertsa)
                qcell_verts.append(ptsa)
                quad_ccents.append(cc1)
                quad_rcircs.append(rc1)
                quad_sa.append(a1)
                quadcell_i.append(ai)
                quad_cents.append(c1)

                self.free_to_merge.remove(ai)

        quad_cells = np.asarray(quad_cells)
        quad_ccents = np.asarray(quad_ccents)
        quad_rcircs = np.asarray(quad_rcircs)
        quad_sa = np.asarray(quad_sa)

        qcell_verts = np.asarray(qcell_verts)
        quad_cents = np.asarray(quad_cents)


        # Reassign all relevant quantities from original tri-mesh:
        self.tri_cells = quad_cells # indices of cells
        self.tri_ccents = quad_ccents # circumcenters
        self.tri_rcircs = quad_rcircs
        self.tri_cents = quad_cents # centroids
        self.tcell_verts = qcell_verts # x,y coordinates of vertices of cells
        self.tri_sa = quad_sa  # surface area of triangle

        self.n_tcell = len(self.tri_cells)  # number of simplexes in trimesh

        self.tri_cell_i = np.asarray([i for i in range(self.n_tcell)])

    def process_primary_edges(self, mesh_version = 'tri'):
        """
        Processes the edges and boundary of the primary mesh.
        :return:
        """

        # Calculate edges for cells trimesh, and the hull:
        all_edges = set()  # set of all vertex pairs
        unique_edges = set()  # set of unique vertex pairs
        hull_points = []
        hull_edges = []
        tri_mids = [] # midpoint of tri edges
        tri_edge_len = [] # length of tri edge
        tri_tang = [] # tri edge tangent

        # Begin by creating a master set of all edges (including duplicate (vi, vj) and (vj, vi) combos:
        for verti_o in self.tri_cells:

            verti_i = np.roll(verti_o, -1)

            for vi, vj in zip(verti_o, verti_i):
                all_edges.add((vi, vj))

        for va, vb in all_edges:

            if (va, vb) in all_edges and (vb, va) not in all_edges:
                # if there isn't a double-pair, then add these edges to the hull:
                # (this is based on the logic that when traversing the points of the
                # triangular simplices, only the boundary edges are traversed once,
                # since they don't have a neighbouring simplex at the bounds.)
                hull_points.append(va)
                hull_points.append(vb)

                hull_edges.append([va, vb])


            if (vb, va) not in unique_edges: # otherwise add the edge to the set
                unique_edges.add((va, vb))


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

        self.bflags_tedges = np.asarray(bflags_tedges)  # indices of edges on the boundary

        self.tri_edge_i = np.linspace(0, self.n_tedges - 1, self.n_tedges, dtype=np.int)


        # Finally, go through and calculate mids, len, and tangents of tri_edges, and prepare a mapping between
        # each vertices and edges:
        for ei, (vi, vj) in enumerate(self.tri_edges):

            # get coordinates associated with each edge
            tpi = self.tri_verts[vi]
            tpj = self.tri_verts[vj]

            tan_t = tpj - tpi

            tri_len = np.linalg.norm(tan_t)

            tan_ti = tan_t/tri_len

            assert (np.round(tri_len, 15) != 0.0), "Tri-edge length equal to zero! Duplicate seed points exist!"

            pptm = (tpi + tpj) / 2  # calculate midpoint of tri-edge:

            tri_mids.append(pptm)
            tri_edge_len.append(tri_len)
            tri_tang.append(tan_ti)

        self.tri_mids = np.asarray(tri_mids)
        self.tri_edge_len = np.asarray(tri_edge_len)
        self.tri_tang = np.asarray(tri_tang)

        self.inner_vvert_i = np.linspace(0, len(self.tri_ccents), len(self.tri_ccents), dtype=np.int)

    def create_tri_map(self):
        """
        Creates the basic mapping between triverts and simplices
        :return:
        """
        # create an array giving a list of simplex indices for each tri_vert
        verts_to_simps = [[] for i in range(len(self.tri_verts))]

        # create a set containing all tri_cell inds as tuples (used to control quad merging):
        self.free_to_merge = set()

        for ci, vertsi_o in enumerate(self.tri_cells):

            self.free_to_merge.add(ci)

            for vi in vertsi_o:
                verts_to_simps[vi].append(ci)

        self.tverts_to_tcell = np.asarray(verts_to_simps) # for each tvert, what simplices does it belong to?

    def create_mappings(self):


        # make the face-to-edge indices mapping for the tri_mesh:
        tri_edges = self.tri_edges.tolist()

        face_to_edges = [[] for ii in range(self.n_tcell)]
        bflags_tcells = []

        # mapping giving list of edge inds for each tri vert
        verts_to_edges = [[] for i in range(len(self.tri_verts))]

        for ci, vertsi_o in enumerate(self.tri_cells):

            vertsi_i = np.roll(vertsi_o, -1)

            for vi, vj in zip(vertsi_o, vertsi_i):

                if [vi, vj] not in tri_edges:
                    # get the index of the opposite sign edge:
                    ea = tri_edges.index([vj, vi])
                    face_to_edges[ci].append(ea) # append the edge index to the array at the cell index
                    verts_to_edges[vi].append(ea) # append the edge index to each vertex ind as well
                    verts_to_edges[vj].append(ea)

                else:
                    # get the forward sign edge:
                    ea = tri_edges.index([vi, vj])
                    face_to_edges[ci].append(ea)
                    verts_to_edges[vi].append(ea) # append the edge index to each vertex ind as well
                    verts_to_edges[vj].append(ea)

                # if any edge is on the boundary, mark the cell
                if ea in self.bflags_tedges:
                    bflags_tcells.append(ci)


        self.tcell_to_tedges = np.asarray(face_to_edges) # tri_face index to tri_edges indices mapping
        self.bflags_tcells = np.asarray(np.unique(bflags_tcells))  # trimesh faces on boundary
        self.tverts_to_tedges = np.asarray(verts_to_edges) # for each tever, what edges does it belong to?

    def map_tvert_to_vorcell(self):
        """
        Create a rudimentary mapping between a trivert and the
        vor cell verts associated with it.

        """
        # For each trivert, find the vor verts for cell surrounding it:
        triverts_to_vorverts = []

        # Points search tree with respect to tri_mids points of edges:
        vor_tree = cKDTree(self.vor_verts)

        # For the set of edges associated with each trivert:
        for vi, edge_inds_o in enumerate(self.tverts_to_tedges):
            edge_inds = np.unique(edge_inds_o) # make the set of edge verts unique!
            mid_pts = self.tri_mids[edge_inds]
            dd, vort_indi_o = vor_tree.query(mid_pts, k=2) # this finds 2 vor 'ridge points' per tri edge
            vort_indis = np.unique(vort_indi_o) # get the unique vor vert inds from set of all ridge inds
            triverts_to_vorverts.append(vort_indis)  # append it to the list of vor verts at the tri vert

        self.triverts_to_vorverts = np.asarray(triverts_to_vorverts)

    def define_vorverts(self):
        """
        When constructing a voronoi diagram, the vor mesh edges intersecting the incredible hull of the
        tri mesh must be defined by a new set of vor mesh vertices added outside the tri mesh hull. This
        method adds these perpendicular bisectors into the vor mesh.

        """

        # Make a search tree of interior vor vertices
        if self.use_centroids:
            vor_tree = cKDTree(self.tri_cents)
            vor_points = self.tri_cents

        else:
            vor_tree = cKDTree(self.tri_ccents)
            vor_points = self.tri_ccents

        vor_verts_bound = [] # vertices at the boundary

        for si, edge_inds in enumerate(self.tcell_to_tedges):

            for ei in edge_inds:

                if ei in self.bflags_tedges:  # If this edge is on the hull/boundary, need to calc bound vorcell points

                    tri_edge_mid = self.tri_mids[ei] # midpoint of tri_edge
                    tri_edge_tang = self.tri_tang[ei] # tangent of tri_edge

                    if self.use_centroids:
                        tri_cell_cent = self.tri_cents[si] # get the ccent or cent point of this tri cell
                    else:
                        tri_cell_cent = self.tri_ccents[si]

                    d_cent_to_mid = np.linalg.norm(tri_edge_mid - tri_cell_cent) # get the distance from cell center
                                                                                 # to the edge midpoint (to use as
                                                                                 # potentially anisotropic radius)

                    # query the interior vor points for the vor index of the nearest neigh vor point to the edge mid:
                    d_va, i_va = vor_tree.query(tri_edge_mid)

                    # get coordinates of the first point of the missing vor ridge at the boundary:
                    vor_pt_a = vor_points[i_va]

                    # get the vor edge tangent vector as the 90o rotation of tri_edge tangent (x --> y, y --> -x):
                    vor_tang = tri_edge_tang*1
                    vor_tang_x = vor_tang[1]
                    vor_tang_y = -vor_tang[0]

                    # calculate the boundary vor verts from hull edge properties:
                    vor_pt_b = vor_pt_a*1

                    vor_pt_b[0] = vor_tang_x*d_cent_to_mid*2 + vor_pt_a[0]
                    vor_pt_b[1] = vor_tang_y*d_cent_to_mid*2 + vor_pt_a[1]

                    vor_verts_bound.append(vor_pt_b)

        # define vor_verts:
        self.vor_verts_bound = np.asarray(vor_verts_bound)

        if self.use_centroids:
            self.vor_verts = np.vstack((self.tri_cents, self.vor_verts_bound))
        else:
            self.vor_verts = np.vstack((self.tri_ccents, self.vor_verts_bound))

    def process_vorcells(self):

        vor_cells = []  # vor_vert inds specifying dual cell for each tri_vert
        vcell_verts = []
        vor_cents = []  # centroid of voronoi polygons
        vor_sa = []  # surface area of voronoi polygons

        tvert_to_vorcell = [] # updated mapping between tverts and vorcells

        for tvi, tedge_inds_o in enumerate(self.tverts_to_tedges):
            # make the edge inds for this trivert unique:
            tedge_inds = np.unique(tedge_inds_o)

            # find the vor verts that match these edges
            vort_indis = np.unique(self.vor_edges[tedge_inds].ravel())

            if len(vort_indis) > 2: # if there's more than 2 vor verts in the dual cell (enough to make a polygon:
                # get voronoi verts corresponding to each voronoi dual-cell simplex:
                vor_pts = self.vor_verts[vort_indis] # unsorted!

                # sort the points counter-clockwise:
                cent = vor_pts.mean(axis=0)  # calculate the centre point
                angles = np.arctan2(vor_pts[:, 1] - cent[1], vor_pts[:, 0] - cent[0])  # calculate point angles
                sorted_region = vort_indis[np.argsort(angles)]  # sort indices counter-clockwise

                vor_cells.append(sorted_region)  # add sorted list to the regions structure

                vpts = self.vor_verts[sorted_region]  # x,y coordinates, sorted
                vcell_verts.append(vpts)

                # Calculate centroid and area of the voronoi polygon:
                cx, cy = self.poly_centroid(vpts)
                vor_cents.append([cx, cy])

                sa = self.area(vpts)
                vor_sa.append(sa)

                # update the map between tvert and vorcell:
                tvert_to_vorcell.append(sorted_region)

            else:
                logs.log_warning("Warning! Non-polygonal Voronoi cell detected!")


        self.triverts_to_vorverts = np.asarray(tvert_to_vorcell)
        self.vor_cells = np.asarray(vor_cells)
        self.vcell_verts = np.asarray(vcell_verts)
        self.vor_cents = np.asarray(vor_cents)
        self.vor_sa = np.asarray(vor_sa)

        self.n_vverts = len(self.vor_verts)
        self.n_vcells = len(self.vor_cells)

        # Define indices arrays for vor mesh:
        self.vor_vert_i = np.linspace(0, self.n_vverts - 1, self.n_vverts, dtype=np.int)
        self.vor_cell_i = np.linspace(0, self.n_vcells - 1, self.n_vcells, dtype=np.int)

        # get final bounds for the cluster:
        xmin = self.vor_verts[:,0].min()
        xmax = self.vor_verts[:,0].max()
        ymin = self.vor_verts[:, 1].min()
        ymax = self.vor_verts[:, 1].max()

        self.xyaxis = [xmin*1.1, xmax*1.1, ymin*1.1, ymax*1.1]

    def process_voredges(self):

        vor_edges = []
        vor_tang = []  # tangent vectors to vor_edges
        vor_edge_len = []  # length of vor_edge
        vor_mids = [] # mids of vor edges

        vor_norm = [] # normals to vor cell surfaces (outwards pointing)
        tri_norm = [] # normals to tri cell surfaces (outwards pointing)
        sflux_n = [] # dot product between voronoi cell surface normal (outwards) and tri-tangent

        # Add in Voronoi points on the outer boundary and define total vor_verts:
        self.define_vorverts()

        # create an initial mapping between triverts and surrounding vor cell verts
        self.map_tvert_to_vorcell()

        # find edges of Voronoi dual mesh:
        # Want vor edges to have the same index as tri_edges and to be
        # perpendicular bisectors; therefore we're going to have one vor_edge vert pair for each tri-edge
        for tei, (ti, tj) in enumerate(self.tri_edges):

            # get the voronoi indices (each sorted counter clockwise) for
            # defining voronoi dual cells surrounding each tri vert point:
            vor_reg_i = self.triverts_to_vorverts[ti]
            vor_reg_j = self.triverts_to_vorverts[tj]

            # find the common inds between the two cells,
            # representing the shared edge:
            # print(vor_reg_i, vor_reg_j)
            shared_ij = np.intersect1d(vor_reg_i, vor_reg_j)


            if len(shared_ij) > 2: # We have a multiple-intersecting vor cell issue to deal with...

                # get tangent vector for tri_edge:
                tan_t = self.tri_tang[tei]

                # get vor vert coords of the polygonal intersection:
                poly_pts = self.vor_verts[shared_ij]

                cent = poly_pts.mean(axis=0)  # calculate the centre point
                angles = np.arctan2(poly_pts[:, 1] - cent[1], poly_pts[:, 0] - cent[0])  # calculate point angles
                sorted_reg = shared_ij[np.argsort(angles)]  # sort indices counter-clockwise

                # get edge pairs for newly defined polygon:
                sorted_reg_i = np.roll(sorted_reg, -1)

                dot_check_vals = [] # array to store dot product checks between each vor edge and tri tangent

                for i, ptsi in enumerate(zip(sorted_reg, sorted_reg_i)):
                    ptsa = self.vor_verts[ptsi[0]]
                    ptsb = self.vor_verts[ptsi[1]]

                    tan_vo = ptsb - ptsa
                    len_v = np.linalg.norm(tan_vo)
                    tan_v = tan_vo / len_v
                    dot_test = np.dot(tan_t, tan_v)

                    dot_check_vals.append(dot_test)

                sorted_dot_check = np.argsort(np.abs(dot_check_vals))

                # redefine shared ij pair in terms of most orthogonal vor edge to tri edge:
                shared_ij = np.asarray([sorted_reg[sorted_dot_check][0], sorted_reg_i[sorted_dot_check][0]])

                vor_edge_ij = shared_ij*1

            elif len(shared_ij) == 2:

                vor_edge_ij = shared_ij*1

            else:
                raise Exception("Empty Voronoi edge.")


            # find points representing vor and tri edges:
            tan_t = self.tri_tang[tei]

            tri_norm.append([tan_t[1], -tan_t[0]])

            vpi = self.vor_verts[vor_edge_ij[0]]
            vpj = self.vor_verts[vor_edge_ij[1]]

            tan_vo = vpj - vpi

            vor_len = np.linalg.norm(tan_vo)

            tan_v = tan_vo/vor_len

            ppvm = (vpi + vpj)/2 # midpoint of voronoi edge
            vor_mids.append(ppvm)

            assert (vor_len != 0.0), "Vor-edge length equal to zero!"

            vor_edge_len.append(vor_len)

            cross_tp = tan_t[0]*tan_v[1] - tan_v[0]*tan_t[1]

            if np.sign(cross_tp) == -1.0:

                vor_edges.append([vor_edge_ij[0], vor_edge_ij[1]])
                vor_tang.append(tan_v)
                norm_v = [tan_v[1], -tan_v[0]]
                vor_norm.append(norm_v)
                sflux_n.append(-(norm_v[0]*tan_t[0] + norm_v[1]*tan_t[1]))

            elif np.sign(cross_tp) == 1.0:

                vor_edges.append([vor_edge_ij[1], vor_edge_ij[0]])
                vor_tang.append(-tan_v)
                norm_v = [-tan_v[1], tan_v[0]]
                vor_norm.append(norm_v)
                sflux_n.append(-(norm_v[0]*tan_t[0] + norm_v[1]*tan_t[1]))

        self.vor_edges = np.asarray(vor_edges)
        self.vor_tang = np.asarray(vor_tang)

        self.vor_edge_len = np.asarray(vor_edge_len)  # length of vor_edge
        self.vor_mids = np.asarray(vor_mids) # edge midpoints

        self.n_vedges = len(self.vor_edges)
        self.vor_edge_i = np.linspace(0, self.n_vedges - 1, self.n_vedges, dtype=np.int)

        self.vor_norm = np.asarray(vor_norm) # normals to vor cell surfaces (outwards pointing)
        self.tri_norm = np.asarray(tri_norm) # normals to tri cell surfaces (outwards pointing)
        self.sflux_n = np.asarray(sflux_n) # dot product between voronoi cell surface normal (outwards) and tri-tangent

    def sanity_check(self):

        # Check to see if some vertices are not used in any simplex:
        unused_tverts = []
        for tvi, sverts in enumerate(self.tverts_to_tcell):
            if len(sverts) == 0:
                unused_tverts.append(tvi)

        unused_tverts = np.asarray(unused_tverts)

        if len(unused_tverts) and self.removed_bad_verts is False:

            logs.log_info("Some seed verts are not used in any mesh simplexes! \n"
                             "Decrease alpha_shape value to avoid this case. \n"
                             "Removing unused seed points from mesh...")

            # if there's a length to the unused tverts, deleted flags from data-structures:
            self.tri_verts = np.delete(self.tri_verts, unused_tverts, axis = 0)
            self.n_tverts = len(self.tri_verts) # number of tri_verts
            self.tri_vert_i = np.linspace(0, self.n_tverts - 1,
                                          self.n_tverts, dtype=np.int)

            # remake the trimesh with reduced set of triverts:
            self.trimesh_core_calcs()

            self.removed_bad_verts = True

        elif len(unused_tverts) and self.removed_bad_verts is True:
            logs.log_warning("Something is seriously wrong with your mesh!")

    def update_metrics(self):
        """
        Updates Voronoi cell area, and tri and vor edge lengths.
        This is required for case of deforming mesh, when
        connectivity of the meshes are not changing.

        """
        # Update tri_ccents, tri_cents as well!
        vor_cents = []  # centroid of voronoi polygons
        vor_sa = []  # surface area of voronoi polygons
        vor_edge_mids = []  # mids of tri cell edges
        vor_edge_len = []  # length of vor_edge
        vor_tang = []  # tangent vectors to vor_edges
        vor_norm = [] # normal vectors to vor_edges


        tri_cents = [] # center of triangular cells #FIXME implement this!
        tri_ccents = [] # circumcenter of triangles
        tri_rcircs = [] # circumradius of triangles
        tri_sa = [] # surface area of triangles
        tri_edge_mids = []  # mids of tri cell edges
        tri_edge_len = []  # length of tri_edge
        tri_tang = []  # tangent vectors to tri_edges
        tri_norm = []  # tangent vectors to vor_edges

        sflux_n = [] # dot product between vor cell normal and tricell tangents


        for pts in self.vcell_verts: # FIXME calculate tri cell properties here too (each vcell index maps to tri vert index)
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

            tan_to = tpj - tpi

            tri_len = np.linalg.norm(tan_to)
            tan_t = tan_to/tri_len

            pptm = (tpi + tpj) / 2  # calculate midpoint of tri-edge:

            tri_edge_mids.append(pptm)

            assert (tri_len != 0.0), "Tri-edge length equal to zero!"

            tri_edge_len.append(tri_len * 1)
            tri_tang.append(1 * tan_t)
            tri_norm.append([tan_t[1], -tan_t[0]])

            vi, vj = self.vor_edges[tei]

            # get x-y points:
            vpi = self.vor_verts[vi]
            vpj = self.vor_verts[vj]

            tan_vo = vpj - vpi

            vor_len = np.linalg.norm(tan_vo)

            tan_v = tan_vo/vor_len

            ppvm = (vpi + vpj)/2 # midpoint of voronoi edge
            vor_edge_mids.append(ppvm)

            assert (vor_len != 0.0), "Vor-edge length equal to zero!"

            vor_edge_len.append(vor_len * 1)
            vor_tang.append(tan_v)
            norm_v = [tan_v[1], -tan_v[0]]
            vor_norm.append(norm_v)

            sflux_n.append(-(norm_v[0] * tan_t[0] + norm_v[1] * tan_t[1]))

        for si, verti in enumerate(self.tri_cells):
            tripts = self.tri_verts[verti]

            cx, cy = self.poly_centroid(tripts)
            tri_cents.append([cx, cy])

            sa = self.area(tripts)  # surface area of triangle
            tri_sa.append(sa)

            if len(verti) == 3:
                vx, vy, r_circ, r_in = self.circumc(tripts[0], tripts[1], tripts[2])
                tri_ccents.append([vx, vy])
                tri_rcircs.append(r_circ)

            elif len(verti) == 4:
                R, area, cx, cy = self.quad_circumc(tripts[0], tripts[1], tripts[2], tripts[3])
                tri_ccents.append([cx, cy])
                tri_rcircs.append(R)

        self.tri_tang = np.asarray(tri_tang)
        self.vor_tang = np.asarray(vor_tang)
        self.tri_cents = np.asarray(tri_cents)
        self.tri_ccents = np.asarray(tri_ccents)
        self.tri_sa = np.asarray(tri_sa)

        self.vor_edge_len = np.asarray(vor_edge_len)  # length of vor_edge
        self.tri_edge_len = np.asarray(tri_edge_len)  # length of tri_edge

        self.tri_mids = np.asarray(tri_edge_mids)  # edge midpoints
        self.vor_mids = np.asarray(vor_edge_mids) # edge midpoints

        self.vor_norm = np.asarray(vor_norm) # normals to vor cell surfaces (outwards pointing)
        self.tri_norm = np.asarray(tri_norm) # normals to tri cell surfaces (outwards pointing)
        self.sflux_n = np.asarray(sflux_n) # dot product between voronoi cell surface normal (outwards) and tri-tang

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

    def create_aux_operators(self):
        """
        Creates auxiliary operators required for curl, vector laplacians, etc. Note these are
        needed for the main mesh, but not for the 'mu-mesh' that is required for tensor work...

        """

        # exterior derivative operator for tri mesh operating on edges to return faces:
        delta_tri_1 = np.zeros((self.n_tcell, self.n_tedges))

        tri_edges = self.tri_edges.tolist()

        for ic, vertis_o in enumerate(self.tri_cells):

            vertis_1 = np.roll(vertis_o, -1)

            for vi, vj in zip(vertis_o, vertis_1):

                if [vi, vj] not in tri_edges:
                    # get the index of the opposite sign edge:
                    ea = tri_edges.index([vj, vi])
                    delta_tri_1[ic, ea] = -1
                else:
                    # get the forward sign edge:
                    ea = tri_edges.index([vi, vj])
                    delta_tri_1[ic, ea] = 1

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

        # ## Note the following creates the delta_vor_1 exterior derivative, however, delta_vor_1 is
        # # not required as: delta_vor_1 = delta_tri_1.T. Keeping only for testing purposes.
        # # exterior derivative operator for vor mesh operating on edges to return faces:
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
        #
        # self.delta_vor_1_inv = np.linalg.pinv(self.delta_vor_1)

        # Create mapping from tri verts to tri centers (Uses Barycentric coordinates to interpolate
        # from verts to circumcentre):

        # M_verts_to_cents = np.zeros((self.n_tcell, self.n_tverts))
        #
        # for ii, edge_inds in enumerate(self.tcell_to_tedges):
        #     a, b, c = self.tri_edge_len[edge_inds] # FIXME don't know that they're triangles!
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
        #     vi, vj, vk = self.tri_cells[ii] # FIXME! Don't know that they're triangles!
        #
        #     M_verts_to_cents[ii, vi] = b1
        #     M_verts_to_cents[ii, vj] = b2
        #     M_verts_to_cents[ii, vk] = b3
        #
        # self.M_verts_to_cents = np.asarray(M_verts_to_cents)
        # self.M_verts_to_cents_inv = np.linalg.pinv(self.M_verts_to_cents)


    #----Mathematical operator functions-----------

    def gradient_xy(self, Sv, gtype = 'tri'):
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

            assert(len(Sv) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            Sd = self.verts_to_verts(Sv, gtype = 'tri') # interpolate to verts of dual mesh

            gS_tri = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, Sv) # gradient with respect to tri mesh
            gS_vor = (1/self.vor_edge_len)*np.dot(self.delta_vor_0, Sd) # gradient with respect to vor mesh

            # gradSx = (1/self.tri_edge_len)*gS*self.tri_tang[:,0]
            # gradSy = (1 / self.tri_edge_len)*gS*self.tri_tang[:, 1]

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            assert(len(Sv) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            Sd = self.verts_to_verts(Sv, gtype = 'vor') # interpolate to verts of dual mesh

            gS_tri = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, Sd) # gradient with respect to tri mesh
            gS_vor = (1/self.vor_edge_len)*np.dot(self.delta_vor_0, Sv) # gradient with respect to vor mesh

            # gradSx = (1/self.vor_edge_len)*gS*self.vor_tang[:,0]
            # gradSy = (1/self.vor_edge_len)*gS*self.vor_tang[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        gradSx = self.tri_tang[:,0]*gS_tri + self.tri_tang[:,1]*gS_vor
        gradSy = self.vor_tang[:,0]*gS_tri + self.vor_tang[:,1]*gS_vor


        return gradSx, gradSy

    def gradient_uv(self, Sv, gtype = 'tri'):
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

            assert(len(Sv) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            gS = np.dot(self.delta_tri_0, Sv) # gradient with respect to tri mesh

            gradSx = (1/self.tri_edge_len)*gS*self.tri_tang[:,0]
            gradSy = (1 / self.tri_edge_len)*gS*self.tri_tang[:, 1]

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            assert(len(Sv) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            gS = np.dot(self.delta_vor_0, Sv) # gradient with respect to vor mesh

            gradSx = (1/self.vor_edge_len)*gS*self.vor_tang[:,0]
            gradSy = (1/self.vor_edge_len)*gS*self.vor_tang[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")


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
            raise Exception("valid gtype is 'tri' or 'vor'")

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
            raise Exception("valid gtype is 'tri' or 'vor'")


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

            sflux_n =  self.sflux_n

            divF = (1/self.vor_sa)*np.dot(-self.delta_tri_0.T, sflux_n*self.vor_edge_len*Ft)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor div!"

            divF = (1/self.tri_sa)*np.dot(self.delta_tri_1, self.tri_edge_len*Ft)

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")


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

            sflux_n = self.sflux_n
            # ensure passed array is of the correct length:
            assert (len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate gradient of S:
            gS = self.gradient(S, gtype='tri')

            # calculate the divergence of the gradient, which is the laplacian:
            lapS = self.div(sflux_n*gS, gtype = 'tri')

        elif gtype == 'vor':

            assert(len(S) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate gradient of S:
            gS = self.gradient(S, gtype='vor')

            # calculate the divergence of the gradient, which is the Laplacian:
            lapS = self.div(gS, gtype = 'vor')

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return lapS

    def lap_inv(self, S, gtype = 'tri'):

        """
        Computes an inverse scalar Laplacian on a scalar variable S as the inverse div of the inverse gradient of the
        scalar.

        If gtype = 'tri', the gradient is taken with respect to the tri_mesh edges with vor_mesh control volumes,
        and if gtype = 'vor', the gradient is taken with respect to the vor_mesh edges with tri_mesh control volumes.

        Note that due to the structure of the grids, the gtype tri lap is closed boundary (zero flux) while the
        gtype vor lap is a 'free'/'open' boundary.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if laplacian is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        lapS  -- the Laplacian of S with 'natural' boundary conditions.

        """
        if gtype == 'tri':

            sflux_n =  self.sflux_n

            # ensure passed array is of the correct length:
            assert(len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate the divergence of the gradient, which is the laplacian:
            lapS_inv = np.dot(self.delta_tri_0_inv,
                              (self.tri_edge_len/
                               (self.vor_edge_len*sflux_n))*np.dot(-self.delta_tri_0_inv.T, S*(self.vor_sa)))

        elif gtype == 'vor':

            # ensure passed array is of the correct length:
            assert(len(S) == len(self.tri_ccents)), "Length of array passed to gradient is not tri_faces length"

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate inverse Laplacian of S:
            lapS_inv = np.dot(self.delta_vor_0_inv,
                   (self.vor_edge_len/self.tri_edge_len)*np.dot(self.delta_tri_1_inv, S*(self.tri_sa)))

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return lapS_inv

    def curl_z(self, Fz, gtype = 'tri'):

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            assert (len(Fz) == self.n_tverts), "Length of array passed to curl is not tri verts length!"


            curlFt = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, Fz)

            # the x, y components of curl are the skew gradient:

            curlFz_x = curlFt*self.vor_tang[:,0]
            curlFz_y = curlFt*self.vor_tang[:,1]


        elif gtype == 'vor':

            assert (len(Fz) == self.n_vverts), "Length of array passed to curl is not vor verts length!"

            curlFt = -(1/self.vor_edge_len)*np.dot(self.delta_vor_0, Fz)

            curlFz_x = curlFt*self.tri_tang[:,0]
            curlFz_y = curlFt*self.tri_tang[:,1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return curlFz_x, curlFz_y

    def curl_z2(self, Fz, gtype = 'tri'):

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            assert (len(Fz) == self.n_tverts), "Length of array passed to curl is not tri verts length!"

            gfx, gfy = self.gradient_xy(Fz, gtype='tri')

            curlFz_x = -gfy
            curlFz_y = gfx


        elif gtype == 'vor':

            assert (len(Fz) == self.n_vverts), "Length of array passed to curl is not vor verts length!"


            gfx, gfy = self.gradient_xy(Fz, gtype='vor')

            curlFz_x = -gfy
            curlFz_y = gfx

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return curlFz_x, curlFz_y

    def curl_xy(self, Fx, Fy, gtype = 'tri'):
        """
        Calculates the curl of an Fx, Fy vector field and returns the curl = Phi_z component.
        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """

        # as there are an equal number of edges in tri and vor grids, check once:
        assert (len(Fx) == self.n_tedges), "Length of array passed to curl_xy is not edges length!"
        assert (len(Fy) == self.n_tedges), "Length of array passed to curl_xy is not edges length!"

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # get tangential component of Fx, Fy with respect to the tri_tangents:
            Ft = Fx*self.tri_tang[:,0] + Fy*self.tri_tang[:,1]

            # calculate the curl (which is a vector in the z-direction with + representing out of page):
            curl_F = (1/self.tri_sa)*np.dot(self.delta_tri_1, (self.tri_edge_len)*Ft)


        elif gtype == 'vor':

            # get tangential component of Fx, Fy with respect to the vor_tangents:
            Ft = Fx*self.vor_tang[:,0] + Fy*self.vor_tang[:,1]

            # calculate the curl (which is a vector in the z-direction with + representing out of page):
            curl_F = (1/self.vor_sa)*np.dot(-self.delta_tri_0.T, (self.vor_edge_len)*Ft)

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return curl_F

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
            raise Exception("valid gtype is 'tri' or 'vor'")

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
            MM_inv = np.abs(self.delta_tri_0.T)*(1/2)

            Sv = np.dot(MM_inv, Sm)

        elif gtype == 'vor':

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"
            assert(len(Sm) == self.n_vedges), "Length of array passed to gradient is not edges length"

            MM_inv = np.abs(self.delta_vor_0.T)*(1/2)

            Sv = np.dot(MM_inv, Sm)

        else:

            raise Exception("valid gtype is 'tri' or 'vor'")

        return Sv

    def verts_to_verts(self, Sv, gtype ='tri', bval = 0.0):
        """
        Maps quantity from vertices of one grid (e.g. tri) to the
        dual (e.g. vor).

        Parameters
        ------------
        Sv    Variable defined on tri or vor mesh verts
        bval  Values for the boundary (converting from tri to vor grid only)
        gtype Type of starting mesh

        Return
        --------
        Sd     Variable defined on opposite mesh to what was input

        """

        if gtype == 'vor':
            Sv_edges = self.verts_to_mids(Sv, gtype='vor')
            path_len_tri = np.dot(np.abs(self.delta_tri_0.T), self.vor_edge_len)
            Sd = np.dot(np.abs(self.delta_tri_0.T), self.vor_edge_len * Sv_edges) / path_len_tri

        elif gtype == 'tri':
            Sv_edges = self.verts_to_mids(Sv, gtype='tri')
            path_len_vor = np.dot(np.abs(self.delta_tri_1), self.tri_edge_len)
            Sd = np.ones(self.n_vverts)*bval
            Sd[self.inner_vvert_i] = np.dot(np.abs(self.delta_tri_1), self.tri_edge_len*Sv_edges)/path_len_vor

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")


        return Sd

    def mids_to_verts2(self, Sm, gtype ='tri', bval = 0.0):
        """
        Maps quantity from mids to vertices of a grid (e.g. tri).

        Parameters
        ------------
        Sv    Variable defined on tri or vor mesh verts
        bval  Values for the boundary (converting from tri to vor grid only)
        gtype Type of starting mesh

        Return
        --------
        Sd     Variable defined on opposite mesh to what was input

        """

        if gtype == 'tri':
            path_len_tri = np.dot(np.abs(self.delta_tri_0.T), self.vor_edge_len)
            Sv = np.dot(np.abs(self.delta_tri_0.T), self.vor_edge_len * Sm) / path_len_tri

        elif gtype == 'vor':
            path_len_vor = np.dot(np.abs(self.delta_tri_1), self.tri_edge_len)
            Sv = np.ones(self.n_vverts)*bval
            Sv[self.inner_vvert_i] = np.dot(np.abs(self.delta_tri_1), self.tri_edge_len*Sm)/path_len_vor

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")


        return Sv

    def curl_of_curl(self, Fz, gtype = 'tri'):
        """
        Calculates the Vector Laplacian for the curl of the curl of a vector field Fx, Fy minus the
        gradient of the divergence of Fx, Fy.

        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # For the Curl of curl, the initial Fz is defined on the opposite points to the mesh you are
            # working with:
            assert(len(Fz) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate gradient of Fz with respect to the vor mesh; the skew gradient is the curl:
            gS = self.gradient(-Fz, gtype='vor')

            # calculate the divergence of the gradient on the vor mesh, which is the curl of the curl on
            # the tri mesh:
            ccS = self.div(gS, gtype = 'vor')

        elif gtype == 'vor':

            sflux_n = self.sflux_n
            # ensure passed array is of the correct length:
            assert (len(Fz) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate gradient of Fz on the tri mesh:
            gS = self.gradient(Fz, gtype='tri')

            # calculate the divergence of the gradient on the tri mesh, which is the curl of the curl on
            # the vor mesh:
            ccS = self.div(sflux_n * gS, gtype='tri')

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return ccS

    def curl_of_curl_inv(self, Fz, gtype = 'tri'):
        """
        Computes an inverse vector Laplacian as the inverse curl of the curl of the z-hat component of a
        vector field Fz.  Note the use of opposite mesh versions (e.g. 'vor' operators when a 'tri' mesh
        is specified, and vice-versa, is *not* a mistake, and is done according to DEC theory.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if laplacian is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        lapS  -- the Laplacian of S with 'natural' boundary conditions.

        """
        if gtype == 'vor':

            sflux_n =  self.sflux_n

            # ensure passed array is of the correct length:
            assert(len(Fz) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate the divergence of the gradient, which is the laplacian:
            ccS_inv = np.dot(self.delta_tri_0_inv,
                              (self.tri_edge_len/
                               (self.vor_edge_len*sflux_n))*np.dot(-self.delta_tri_0_inv.T, Fz*(self.vor_sa)))

        elif gtype == 'tri':

            # ensure passed array is of the correct length:
            assert(len(Fz) == len(self.tri_ccents)), "Length of array passed to gradient is not tri_faces length"

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate inverse Laplacian of S:
            ccS_inv = np.dot(self.delta_vor_0_inv,
                   (self.vor_edge_len/self.tri_edge_len)*np.dot(-self.delta_tri_1_inv, Fz*(self.tri_sa)))

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return ccS_inv

    def grad_of_div(self, Fx, Fy, gtype = 'tri'):
        """
        Calculates the gradient of the divergence of a vector field Fx, Fy.

        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """
        # as there are an equal number of edges in tri and vor grids, check once:
        assert (len(Fx) == self.n_tedges), "Length of array passed to gradient is not edges length!"
        assert (len(Fy) == self.n_tedges), "Length of array passed to gradient is not edges length!"

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # get tangential component of Fx, Fy with respect to the tri_tangents:
            Ft = Fx*self.tri_tang[:,0] + Fy*self.tri_tang[:,1]

            div_F = (1/self.vor_sa)*np.dot(-self.delta_tri_0.T, self.vor_edge_len*Ft)

            gdivF = (1/self.tri_edge_len)*np.dot(self.delta_tri_0, div_F)

            gdivFx = gdivF*self.tri_tang[:,0]
            gdivFy = gdivF*self.tri_tang[:, 1]


        elif gtype == 'vor':

            # get tangential component of Fx, Fy with respect to the vor_tangents:
            Ft = Fx*self.vor_tang[:,0] + Fy*self.vor_tang[:,1]


            div_F = (1/self.tri_sa)*np.dot(self.delta_tri_1, self.tri_edge_len*Ft)

            gdivF = (1/self.vor_edge_len)*np.dot(self.delta_vor_0, div_F)

            gdivFx = gdivF*self.vor_tang[:,0]
            gdivFy = gdivF*self.vor_tang[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return gdivFx, gdivFy

    def grad_of_div_inv(self, Fx, Fy, gtype = 'tri'):
        """
        Calculates the inverse gradient of the divergence of a vector field Fx, Fy.

        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """
        # as there are an equal number of edges in tri and vor grids, check once:
        assert (len(Fx) == self.n_tedges), "Length of array passed to gradient is not edges length!"
        assert (len(Fy) == self.n_tedges), "Length of array passed to gradient is not edges length!"

        # Vector Laplacians can only be computed for case where all operators are made:
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # get tangential component of Fx, Fy with respect to the tri_tangents:
            Ft = Fx*self.tri_tang[:,0] + Fy*self.tri_tang[:,1]

            gdivF_inv = (1/self.vor_edge_len)*np.dot(-self.delta_tri_0_inv.T,
                                                     self.vor_sa*np.dot(self.delta_tri_0_inv,
                                                                        self.tri_edge_len*Ft))

            gdivFx = gdivF_inv*self.tri_tang[:,0]
            gdivFy = gdivF_inv*self.tri_tang[:, 1]



        elif gtype == 'vor':

            # get tangential component of Fx, Fy with respect to the vor_tangents:
            Ft = Fx*self.vor_tang[:,0] + Fy*self.vor_tang[:,1]


            gdivF_inv = (1/self.tri_edge_len)*np.dot(self.delta_tri_1_inv,
                                                     self.tri_sa*np.dot(self.delta_vor_0_inv,
                                                                        self.vor_edge_len*Ft))

            gdivFx = gdivF_inv*self.vor_tang[:,0]
            gdivFy = gdivF_inv*self.vor_tang[:,1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return gdivFx, gdivFy

    def curl_of_curl_xy(self, Fx, Fy, gtype = 'tri'):
        """
        Calculates the curl of the curl of a vector field Fx, Fy

        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """
        # as there are an equal number of edges in tri and vor grids, check once:
        assert (len(Fx) == self.n_tedges), "Length of array passed to gradient is not edges length!"
        assert (len(Fy) == self.n_tedges), "Length of array passed to gradient is not edges length!"

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # get tangential component of Fx, Fy with respect to the tri_tangents:
            Ft = Fx*self.tri_tang[:,0] + Fy*self.tri_tang[:,1]

            # calculate the curl of the curl:
            ccft = (1/self.vor_edge_len)*np.dot(self.delta_vor_0,
                                                          (1/self.tri_sa)*np.dot(self.delta_tri_1,
                                                                                 (self.tri_edge_len)*Ft))

            lapFx = ccft*self.tri_tang[:,0]
            lapFy = ccft*self.tri_tang[:, 1]


        elif gtype == 'vor': #FIXME this may be off by a negative sign -- untested!

            # get tangential component of Fx, Fy with respect to the vor_tangents:
            Ft = Fx*self.vor_tang[:,0] + Fy*self.vor_tang[:,1]

            # calculate the curl of the curl:
            ccft = -(1/self.tri_edge_len)*np.dot(self.delta_tri_0,
                                                          (1/self.vor_sa)*np.dot(-self.delta_tri_0.T,
                                                                                 (self.vor_edge_len)*Ft))

            lapFx = ccft*self.vor_tang[:,0]
            lapFy = ccft*self.vor_tang[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return lapFx, lapFy

    def curl_of_curl_xy_inv(self, Fx, Fy, gtype = 'tri'):
        """
        Calculates the inverse Vector Laplacian for a *divergence free* field Fx, Fy as the inverse of the
        curl of the curl of Fx, Fy. Returns components lapFx_inv and lapFy_inv tangential to the mesh gtype used.

        :param Fx:
        :param Fy:
        :param gtype:
        :return:
        """
        # as there are an equal number of edges in tri and vor grids, check once:
        assert (len(Fx) == self.n_tedges), "Length of array passed to gradient is not edges length!"
        assert (len(Fy) == self.n_tedges), "Length of array passed to gradient is not edges length!"

        # Vector Laplacians can only be computed for
        assert (self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

        if gtype == 'tri':

            # get tangential component of Fx, Fy with respect to the tri_tangents:
            Ft = Fx * self.tri_tang[:, 0] + Fy * self.tri_tang[:, 1]

            # calculate the inverse curl of the curl:
            lapFt_inv = (1/self.tri_edge_len)*np.dot(self.delta_tri_1_inv,
                                                     self.tri_sa*np.dot(self.delta_vor_0_inv,
                                                                        Ft*self.vor_edge_len))

            lapFx_inv = lapFt_inv*self.tri_tang[:, 0]
            lapFy_inv = lapFt_inv*self.tri_tang[:, 1]


        elif gtype == 'vor':

            # get tangential component of Fx, Fy with respect to the vor_tangents:
            Ft = Fx * self.vor_tang[:, 0] + Fy * self.vor_tang[:, 1]

            # calculate the inverse curl of the curl:
            lapFt_inv = (1/self.vor_edge_len)*np.dot(-self.delta_tri_0_inv.T,
                                                     self.vor_sa*np.dot(self.delta_tri_0_inv,
                                                                        Ft*self.tri_edge_len))

            lapFx_inv = lapFt_inv*self.vor_tang[:, 0]
            lapFy_inv = lapFt_inv*self.vor_tang[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return lapFx_inv, lapFy_inv

    def helmholtz_hodge(self, Fx, Fy, gtype = 'tri'):
        """
        Decomposes a vector field Fx, Fy into curl-free (gPhi_x, gPhi_y) and div-free (cPsi_x, cPsi_y) components
        using the Helmholtz-Hodge (HH) decomposition, where F is constructed in terms of the gradient of a scalar
        potential Phi and the curl of a (z-axis oriented) vector potential Psi_z:

        F = grad Phi + curl Psi_z + H

        (where 'H' is a harmonic field component that can be found after the HH decomosition is complete by
         calculating: H = F - grad Phi - curl Psi_z.)

        Note that because natural connectivity exists in the orthogonal DEC grids, performing the HH decomp
        with gtype = 'tri' represents an 'open' or 'free' boundary condition for the fields (especially evident
        in the div-free component), whereas gtype = 'vor' represents a closed boundary where there is zero normal
        flux of the div-free component.

        Parameters
        -----------
        Fx: Force field x component (on DEC mesh mids)
        Fy: Force field y component (on DEC mesh mids)
        gtype: type of mesh to work with

        Returns
        ----------
        cPsi_x, cPsi_y: curl of vector potential Psi_z; div-free components of initial field F (on DEC mesh mids)
        gPhi_x, gPhi_y: grad of scalar potential Phi; curl-free components of initial field F (on DEC mesh mids)

        """
        if gtype == 'tri':
            op_gtype = 'vor'

        else:
            op_gtype = 'tri'
        # Solving for the curl-free vector field:
        # take the divergence of the vector field:
        divF = self.div_xy(Fx, Fy, gtype=gtype)

        # The scalar potential of the curl-free component is given by:
        Phi = self.lap_inv(divF, gtype=gtype)

        # Where the curl-free vector field is given by:
        gPhi_x, gPhi_y = self.gradient_uv(Phi, gtype=gtype)

        # Solving for the divergence-free vector field:
        # take the curl of the vector field:
        curlF = self.curl_xy(Fx, Fy, gtype=gtype)

        # The vector potential of the div-free component is given by:
        Psi_z = self.curl_of_curl_inv(curlF, gtype=gtype)

        # Where the div-free component of Fx, Fy is given by:
        cPsi_x, cPsi_y = self.curl_z(Psi_z, gtype=op_gtype)

        return cPsi_x, cPsi_y, gPhi_x, gPhi_y

    def biharmonic(self, S, gtype = 'tri'):
        if gtype == 'tri':

            sflux_n = self.sflux_n
            # ensure passed array is of the correct length:
            assert (len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate gradient of S:
            gS = self.gradient(S, gtype='tri')

            # calculate the divergence of the gradient, which is the laplacian:
            lapS = self.div(sflux_n*gS, gtype = 'tri')

            glapS = self.gradient(lapS, gtype = 'tri')

            bihS = self.div(glapS*self.sflux_n, gtype='tri')

        elif gtype == 'vor':

            assert(len(S) == self.n_vverts), "Length of array passed to gradient is not vor_verts length"

            assert(self.make_all_operators), "This mesh hasn't computed auxillary operators to calculate vor grad"

            # calculate gradient of S:
            gS = self.gradient(S, gtype='vor')

            # calculate the divergence of the gradient, which is the Laplacian:
            lapSo = self.div(gS, gtype = 'vor')

            lapS = np.zeros(len(self.vor_verts))
            lapS[self.inner_vvert_i] = lapSo

            glapS = self.gradient(lapS, gtype='vor')

            bihS = self.div(glapS, gtype='vor')

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return bihS

    def biharmonic_inv(self, S, gtype = 'tri'):

        """
        Computes an inverse scalar Laplacian on a scalar variable S as the inverse div of the inverse gradient of the
        scalar.

        If gtype = 'tri', the gradient is taken with respect to the tri_mesh edges with vor_mesh control volumes,
        and if gtype = 'vor', the gradient is taken with respect to the vor_mesh edges with tri_mesh control volumes.

        Note that due to the structure of the grids, the gtype tri lap is closed boundary (zero flux) while the
        gtype vor lap is a 'free'/'open' boundary.

        Parameters
        -----------
        S   -- a scalar array defined on tri_verts or vor_verts, depending on gtype
        gtype -- specifies if laplacian is taken with respect to tri mesh or vor mesh

        Returns
        ----------
        lapS  -- the Laplacian of S with 'natural' boundary conditions.

        """
        if gtype == 'tri':

            # ensure passed array is of the correct length:
            assert(len(S) == self.n_tverts), "Length of array passed to gradient is not tri_verts length"

            # calculate the divergence of the gradient, which is the laplacian:
            lapinvS = self.lap_inv(S, gtype='tri')

            bihS_inv = self.lap_inv(lapinvS, gtype = 'tri')

        elif gtype == 'vor':

            # ensure passed array is of the correct length:
            assert(len(S) == len(self.tri_ccents)), "Length of array passed to gradient is not tri_faces length"

            # calculate the divergence of the gradient, which is the laplacian:
            lapinvSo = self.lap_inv(S, gtype='vor')

            lapinvS = lapinvSo[self.inner_vvert_i]
            # lapinvS[self.inner_vvert_i] = lapinvSo

            bihS_inv = self.lap_inv(lapinvS, gtype='vor')

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        return bihS_inv

    #----Mesh Refinement-------------------------------

    def calc_tri(self):

        if self.mesh_type == 'tri':
            self.tri_ccents = []
            self.tri_rcircs = [] # circumradius of triangle
            self.tri_sa = []  # surface area of triangle faces

            for i, vert_inds in enumerate(self.tri_cells):

                abc = self.tri_verts[vert_inds]   # FIXME! These may be quad cells!!
                vx, vy, r_circ, r_in = self.circumc(abc[0], abc[1], abc[2])

                sa = self.area(abc)  # surface area of triangle

                self.tri_ccents.append([vx, vy])
                self.tri_rcircs.append(r_circ)  # circumradius of triangle
                self.tri_sa.append(sa)  # surface area of triangle faces

            self.tri_ccents = np.asarray(self.tri_ccents)
            self.tri_rcircs = np.asarray(self.tri_rcircs) # circumradius of triangle
            self.tri_sa = np.asarray(self.tri_sa)  # surface area of triangle faces

    def mesh_quality_calc(self):

        # calculate quality of the mesh in terms of the difference between tri cents and tri ccents

        uu = (np.linalg.norm(self.tri_ccents - self.tri_cents, axis = 1))**2

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

    def refine_mesh(self, max_steps=25, convergence=7.5, fix_bounds = True):

        # if self.mesh_type == 'tri':

        logs.log_info("Initializing Voronoi mesh optimization...")

        opti_steps = np.arange(max_steps)

        ui = self.mesh_quality_calc()

        UU = np.sum(ui)/self.cell_radius**2

        for i in opti_steps:

            self.removed_bad_verts = False  # reset flag for empty tri_vert removal

            if UU > convergence:

                self.removed_bad_verts = False

                # Continuously reassign tri_verts to vor_centres, without affecting the boundary
                if fix_bounds:
                    self.tri_verts[self.biocell_i] = self.vor_cents[self.biocell_i]*1

                # Continuously reassign tri_verts to vor_centres, affecting the boundary
                else:
                    self.tri_verts = self.vor_cents*1

                self.pre_mesh()

                ui = self.mesh_quality_calc()

                UU = np.sum(ui)/self.cell_radius**2

                conv_mess = "Step {}: mesh energy {}".format(i, UU)
                logs.log_info(conv_mess)

            else:
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

        self.removed_bad_verts = False # reset flag for empty tri_vert removal

        for ii, (poly_ind, cell_poly, vor_cent) in enumerate(zip(self.vor_cells,
                                                                 self.vcell_verts,
                                                                 self.vor_cents)):

            # see if voronoi cell center is within the clip curve boundary:
            cent_val = imagemask.clipping_function(vor_cent[0], vor_cent[1])
            if cent_val != 0.0:
                cent_val = 1.0

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

                # the region's points are in the clipping func range, and the voronoi cent lays in the bound:
                elif point_check.sum() > 0.0 and point_check.sum() < len(
                        cell_poly) and cent_val == 1.0:

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

        self.pre_mesh()

    def laplacian_smoothing(self, stepsize):
        """
        Laplacian smoothing of the mesh, using an Implicit Euler
        updating scheme.

        stepsize: How much smoothing to apply (2.0e-6 x cell_radius is good)

        """
        self.removed_bad_verts = False # reset flag for empty tri_vert removal

        II = np.diag(np.ones(self.n_tverts)) # Identity matrix
        HH2 = np.diag(1 / self.vor_sa)  # Hodge star 20
        HH1 = np.diag(self.vor_edge_len / self.tri_edge_len) # Hodge star 11

        term1 = np.dot(HH2, -self.delta_tri_0.T)
        term2 = np.dot(HH1, self.delta_tri_0)

        LL = np.dot(term1, term2) # Forwards Laplacian operator

        MM = (II - stepsize * LL)    # Matrix equation from diffusion equation
        MM_inv = np.linalg.pinv(MM)  # Least-squares solution matrix

        self.tri_verts[:,0] = np.dot(MM_inv, self.tri_verts[:,0]) # Implicit Euler update solution
        self.tri_verts[:,1] = np.dot(MM_inv, self.tri_verts[:,1])


        # # # Laplacian smoothing of the mesh using explicit Euler:
        # delx = self.lap(self.tri_verts[:, 0])
        # dely = self.lap(self.tri_verts[:, 1])
        #
        # self.tri_verts[:, 0] = self.tri_verts[:, 0] + delx*stepsize
        # self.tri_verts[:, 1] = self.tri_verts[:, 1] + dely*stepsize



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

    def quad_circumc(self, A, B, C, D):

        # calculate lengths of all sides
        a = np.linalg.norm(B - A)
        b = np.linalg.norm(C - B)
        c = np.linalg.norm(D - C)
        d = np.linalg.norm(A - D)

        # semiperimeter:
        s = (1 / 2) * (a + b + c + d)

        area = np.sqrt((s - a) * (s - b) * (s - c) * (s - d))

        R = (1 / 4) * (np.sqrt((a * c + b * d) * (a * d + b * c) * (a * b + c * d)) / area)

        # Calculate the diagonal tangent:
        tando = C - A
        tanl = np.linalg.norm(C - A)
        tand = tando / tanl

        # cicumcentre:
        cx = A[0] + (R) * tand[0]
        cy = A[1] + (R) * tand[1]

        return R, area, cx, cy

    #----Tests of DEC computations--------------------

    def plot_test_A(self, a=0.02, b=5.0e-6, gtype = 'vor', size = (10, 8), print_errors = True):
        """
        Generates an analytical test function with analytical gradient, laplacian,
        and curl for comparison with discrete calculations.

        As the 'vor' mesh represents an open/free boundary, this condition best
        replicates the analytical math equations.

        """

        # Generate analytical math:

        if gtype == 'tri':
            xo = self.tri_verts[:, 0]
            yo = self.tri_verts[:, 1]

        elif gtype == 'vor':
            xo = self.vor_verts[:, 0]
            yo = self.vor_verts[:, 1]

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        xm = self.tri_mids[:,0]
        ym = self.tri_mids[:,1]

        # Test function on vor_verts:
        self.foo = np.sin((a * xo * yo) / b ** 2)

        amin = self.foo.min()
        amax = self.foo.max()

        self.gfxm = ((a * ym) / b ** 2) * np.cos((a * xm * ym) / b ** 2)
        self.gfym = ((a * xm) / b ** 2) * np.cos((a * xm * ym) / b ** 2)

        self.gfx = ((a * yo) / b ** 2) * np.cos((a * xo * yo) / b ** 2)
        self.gfy = ((a * xo) / b ** 2) * np.cos((a * xo * yo) / b ** 2)

        self.gf_mag = np.sqrt(self.gfx ** 2 + self.gfy ** 2)
        # self.gf_mmag = self.verts_to_mids(self.gf_mag, gtype=gtype)
        self.gf_mmag = np.sqrt(self.gfxm**2 + self.gfym**2)

        bmin = self.gf_mag.min()
        bmax = self.gf_mag.max()

        # second derivatives
        self.d2fx = -(((a * yo) / b ** 2) ** 2) * np.sin((a * xo * yo) / b ** 2)
        self.d2fy = -(((a * xo) / b ** 2) ** 2) * np.sin((a * xo * yo) / b ** 2)

        self.div_foo = self.d2fx + self.d2fy

        cmin = self.div_foo.min()
        cmax = self.div_foo.max()

        # Curl of foo as phi_z = foo
        self.cfx = self.gfym
        self.cfy = -self.gfxm
        self.cf_mag = np.sqrt(self.cfx** 2 + self.cfy**2)


        # Generate DEC-computed quantities (denoted by '_')-------------------------------------

        # Gradient at edges, in xy components, and div, also using xy components:
        self.gfx_, self.gfy_ = self.gradient_xy(self.foo, gtype=gtype)
        self.div_foo_ = self.div_xy(self.gfx_, self.gfy_, gtype=gtype)

        # # Interpolate gradient to verts and get components
        # self.gfx_ = self.mids_to_verts(gfxmv, gtype=gtype)
        # self.gfy_ = self.mids_to_verts(gfymv, gtype=gtype)
        self.gf_mag_ = np.sqrt(self.gfx_**2 + self.gfy_** 2)

        self.lap_foo_ = self.lap(self.foo, gtype=gtype)

        if gtype == 'vor':
            self.lap_inv_ = self.lap_inv(self.div_foo[self.inner_vvert_i], gtype=gtype)

        else:
            self.div_foo[self.bflags_tverts] = 0.0
            self.lap_inv_ = self.lap_inv(self.div_foo, gtype=gtype)

        # Curl on the vor grid with Phi_z = foov
        self.cfx_, self.cfy_ = self.curl_z2(self.foo, gtype=gtype)
        # self.cfx_ = self.mids_to_verts(cfxm, gtype = gtype)
        # self.cfy_ = self.mids_to_verts(cfym, gtype = gtype)
        self.cf_mag_ = np.sqrt(self.cfx_**2 + self.cfy_**2)

        # RMS errors between quantities:
        self.error_grad = np.sqrt((self.gf_mmag - self.gf_mag_)**2)

        if gtype == 'vor':
            self.error_div = np.sqrt((self.div_foo[self.inner_vvert_i] - self.div_foo_)**2)
            self.error_lap = np.sqrt((self.div_foo[self.inner_vvert_i] - self.lap_foo_)**2)

        else:
            self.error_div = np.sqrt((self.div_foo - self.div_foo_) ** 2)
            self.error_lap = np.sqrt((self.div_foo - self.lap_foo_) ** 2)

        self.error_lap_inv = np.sqrt((self.foo - self.lap_inv_)**2)
        self.error_curl = np.sqrt((self.cf_mag - self.cf_mag_)**2)

        self.test_errors = {'grad': self.error_grad.mean()/np.abs(self.gf_mag).max(),
                            'div': self.error_div.mean()/np.abs(self.div_foo).max(),
                            'lap': self.error_lap.mean()/np.abs(self.div_foo).max(),
                            'lap_inv': self.error_lap_inv.mean()/np.abs(self.foo).max(),
                            'curl': self.error_curl.mean()/np.abs(self.cf_mag).max()}

        if print_errors:
            print("Percent error on estimates:")
            for k, v in self.test_errors.items():
                print(k, np.round(v*100, 4))

        # Plot the results:
        fig, axarr = plt.subplots(2, 3, figsize=size)

        axarr[0, 0].xaxis.set_ticklabels([])
        axarr[0, 0].yaxis.set_ticklabels([])

        axarr[0, 1].xaxis.set_ticklabels([])
        axarr[0, 1].yaxis.set_ticklabels([])

        axarr[0, 2].xaxis.set_ticklabels([])
        axarr[0, 2].yaxis.set_ticklabels([])

        axarr[1, 0].xaxis.set_ticklabels([])
        axarr[1, 0].yaxis.set_ticklabels([])

        axarr[1, 1].xaxis.set_ticklabels([])
        axarr[1, 1].yaxis.set_ticklabels([])

        axarr[1, 2].xaxis.set_ticklabels([])
        axarr[1, 2].yaxis.set_ticklabels([])

        tp1 = axarr[0, 0].tripcolor(xo, yo, self.foo)
        cb1 = fig.colorbar(tp1, ax=axarr[0, 0], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[0, 0].set_title('F (exact)')
        axarr[0, 0].axis('equal')
        axarr[0, 0].axis('off')

        # ---------------------

        tp2 = axarr[0, 1].tripcolor(xm, ym, self.gf_mmag)
        cb2 = fig.colorbar(tp2, ax=axarr[0, 1], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[0, 1].quiver(xo, yo, self.gfx, self.gfy, headwidth = 5)
        axarr[0, 1].set_title(r'$\nabla F (exact)$')
        axarr[0, 1].axis('equal')
        axarr[0, 1].axis('off')

        # ---------------------
        tp3 = axarr[0, 2].tripcolor(xo, yo, self.div_foo)
        cb3 = fig.colorbar(tp3, ax=axarr[0, 2], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[0, 2].set_title(r'$\nabla^{2}F$ (exact)')
        axarr[0, 2].axis('equal')
        axarr[0, 2].axis('off')

        # ---------------------

        tp4 = axarr[1, 0].tripcolor(xo, yo, self.lap_inv_)
        cb4 = fig.colorbar(tp4, ax=axarr[1, 0], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[1, 0].set_title(r'$(\nabla^{2})^{-1}(\nabla\cdot\nabla F)_{exact} (DEC)$')
        axarr[1, 0].axis('equal')
        axarr[1, 0].axis('off')

        tp5 = axarr[1, 1].tripcolor(xm, ym, self.gf_mag_)
        cb5 = fig.colorbar(tp5, ax=axarr[1, 1], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[1, 1].quiver(xm, ym, self.gfx_, self.gfy_, headwidth = 5)
        axarr[1, 1].set_title(r'$\nabla F (DEC)$')
        axarr[1, 1].axis('equal')
        axarr[1, 1].axis('off')

        if gtype == 'vor':
            if self.use_centroids:
                tp6 = axarr[1, 2].tripcolor(self.tri_cents[:, 0], self.tri_cents[:, 1], self.lap_foo_)

            else:
                tp6 = axarr[1, 2].tripcolor(self.tri_ccents[:, 0], self.tri_ccents[:, 1], self.lap_foo_)

        else:
            tp6 = axarr[1, 2].tripcolor(xo, yo, self.lap_foo_)

        cb6 = fig.colorbar(tp6, ax=axarr[1, 2], orientation='horizontal',
                           fraction=0.025, pad=0.01)
        axarr[1, 2].set_title(r'$\nabla^{2}F$ (DEC)')
        axarr[1, 2].axis('equal')
        axarr[1, 2].axis('off')

        tick_locator = ticker.MaxNLocator(nbins=1)

        cb1.locator = tick_locator
        cb1.update_ticks()

        cb2.locator = tick_locator
        cb2.update_ticks()

        cb3.locator = tick_locator
        cb3.update_ticks()

        cb4.locator = tick_locator
        cb4.update_ticks()

        cb5.locator = tick_locator
        cb5.update_ticks()

        cb6.locator = tick_locator
        cb6.update_ticks()

        fig.subplots_adjust(wspace=0.0)

        return fig, axarr

    def plot_test_B(self, a=0.02, b=5.0e-6, size = (10, 6), gtype = 'vor', print_errors = True):
        """
        Compare DEC-based calculation of gradient of gradient with an
        analytical solution.

        Also, all vector fields are ultimately mapped to 'tri' verts because this makes for
        a clearer plot, therefore there is no mesh option for this test.

        """

        # xo = self.tri_verts[:, 0]
        # yo = self.tri_verts[:, 1]

        if gtype == 'tri':
            xo = self.tri_verts[:, 0]
            yo = self.tri_verts[:, 1]
            opgtype = 'vor' # opposite (dual) grid

        elif gtype == 'vor':
            xo = self.vor_verts[:, 0]
            yo = self.vor_verts[:, 1]
            opgtype = 'tri'  # opposite (dual) grid

        else:
            raise Exception("valid gtype is 'tri' or 'vor'")

        xm = self.tri_mids[:, 0]
        ym = self.tri_mids[:, 1]

        # Test function on vor_verts:
        self.foo = np.sin((a * xo * yo) / b ** 2)

        self.gfx = ((a * ym) / b ** 2) * np.cos((a * xm * ym) / b ** 2)
        self.gfy = ((a * xm) / b ** 2) * np.cos((a * xm * ym) / b ** 2)

        self.gf_mag = np.sqrt(self.gfx ** 2 + self.gfy ** 2)

        # Gradient of the gradient:
        self.gfxx = -(((a * ym) / b ** 2) ** 2) * np.sin((a * xm * ym) / b ** 2)
        self.gfxy = (a * np.cos((a * xm * ym) / b ** 2)) / b ** 2 - (
                    a ** 2 * xm * ym * np.sin((a * xm * ym) / b ** 2)) / b ** 4
        self.gfyx = (a * np.cos((a * xm * ym) / b ** 2)) / b ** 2 - (
                    a ** 2 * xm * ym * np.sin((a * xm * ym) / b ** 2)) / b ** 4
        self.gfyy = -(((a * xm) / b ** 2) ** 2) * np.sin((a * xm * ym) / b ** 2)

        #---------
        # Discrete calcs

        # Gradient at edges, in xy components, and div, also using xy components:
        self.gfx_, self.gfy_ = self.gradient_xy(self.foo, gtype=gtype)

        # Interpolate gradient to verts of the dual mesh and get components
        self.gfxi = self.mids_to_verts(self.gfx_, gtype=gtype)
        self.gfyi = self.mids_to_verts(self.gfy_, gtype=gtype)

        # # Interpolate gradient to verts of the original mesh (for plotting)
        # self.gfx_ = self.mids_to_verts(gfxmv, gtype=gtype)
        # self.gfy_ = self.mids_to_verts(gfymv, gtype=gtype)

        # Gradient of the gradient taken from the dual mesh (from verts mapping):
        self.gfxx_, self.gfxy_ = self.gradient_xy(self.gfxi, gtype=gtype)
        self.gfyx_, self.gfyy_ = self.gradient_xy(self.gfyi, gtype=gtype)

        # # each component needs to be mapped to verts (of tri grid type for nice plotting):
        # self.gfxx_ = self.mids_to_verts(gfxxm, gtype=gtype)
        # self.gfxy_ = self.mids_to_verts(gfxym, gtype=gtype)
        # self.gfyx_ = self.mids_to_verts(gfyxm, gtype=gtype)
        # self.gfyy_ = self.mids_to_verts(gfyym, gtype=gtype)

        self.error_gfx = np.sqrt((self.gfx - self.gfx_) ** 2)
        self.error_gfy = np.sqrt((self.gfy - self.gfy_) ** 2)
        self.error_gfxx = np.sqrt((self.gfxx - self.gfxx_)**2)
        self.error_gfxy = np.sqrt((self.gfxy - self.gfxy_)**2)
        self.error_gfyx = np.sqrt((self.gfyx - self.gfyx_)**2)
        self.error_gfyy = np.sqrt((self.gfyy - self.gfyy_)**2)

        self.test_errors = {'gfxx': self.error_gfxx.mean()/np.abs(self.gfxx).max(),
                            'gfxy': self.error_gfxy.mean()/np.abs(self.gfxy).max(),
                            'gfyx': self.error_gfyx.mean()/np.abs(self.gfyx).max(),
                            'gfyy': self.error_gfyy.mean()/np.abs(self.gfyy).max(),
                            'gfx': self.error_gfx.mean() / np.abs(self.gfx).max(),
                            'gfy': self.error_gfy.mean() / np.abs(self.gfy).max(),
                            }

        if print_errors:
            print("Percent error on estimates:")
            for k, v in self.test_errors.items():
                print(k, np.round(v*100, 4))

        # Plot the results:
        fig, axarr = plt.subplots(2, 3, figsize=size)

        axarr[0, 0].xaxis.set_ticklabels([])
        axarr[0, 0].yaxis.set_ticklabels([])

        axarr[0, 1].xaxis.set_ticklabels([])
        axarr[0, 1].yaxis.set_ticklabels([])

        axarr[0, 2].xaxis.set_ticklabels([])
        axarr[0, 2].yaxis.set_ticklabels([])

        axarr[1, 0].xaxis.set_ticklabels([])
        axarr[1, 0].yaxis.set_ticklabels([])

        axarr[1, 1].xaxis.set_ticklabels([])
        axarr[1, 1].yaxis.set_ticklabels([])

        axarr[1, 2].xaxis.set_ticklabels([])
        axarr[1, 2].yaxis.set_ticklabels([])

        axarr[0, 0].quiver(xm, ym, self.gfx, self.gfx, color = 'purple', zorder = 10, headwidth = 5)
        axarr[0, 0].set_title(r'$\nabla\,F\,(exact)$ ')
        # axarr[0, 0].axis('tight')
        axarr[0, 0].axis('off')
        axarr[0, 0].axis('equal')

        axarr[0, 1].quiver(xm, ym, self.gfxx, self.gfxy, color = 'red', zorder = 10, headwidth = 5)
        axarr[0, 1].set_title(r'$\nabla_{x}\nabla_{x}\,Fx (exact)$')
        # axarr[0, 1].axis('tight')
        axarr[0, 1].axis('off')
        axarr[0, 1].axis('equal')

        axarr[0, 2].quiver(xm, ym, self.gfyx, self.gfyy, color = 'blue', zorder = 10, headwidth = 5)
        axarr[0, 2].set_title(r'$\nabla_{y}\nabla_{x}\,F (exact) $')
        # axarr[0, 2].axis('tight')
        axarr[0, 2].axis('off')
        axarr[0, 2].axis('equal')

        axarr[1, 0].quiver(xm, ym, self.gfx_, self.gfy_, color = 'purple', zorder = 10, headwidth = 5)
        axarr[1, 0].set_title(r'$\nabla\,F\,(DEC)$ ')
        # axarr[1, 0].axis('tight')
        axarr[1, 0].axis('off')
        axarr[1, 0].axis('equal')

        axarr[1, 1].quiver(xm, ym, self.gfxx_, self.gfxy_, color = 'red', zorder = 10, headwidth = 5)
        axarr[1, 1].set_title(r'$\nabla_{x}\nabla_{x}\,Fx (DEC)$')
        # axarr[1, 1].axis('tight')
        axarr[1, 1].axis('off')
        axarr[1, 1].axis('equal')

        axarr[1, 2].quiver(xm, ym, self.gfyx_, self.gfyy_, color = 'blue', zorder = 10, headwidth = 5)
        axarr[1, 2].set_title(r'$\nabla_{y}\nabla_{x}\,F (DEC) $')
        # axarr[1, 2].axis('tight')
        axarr[1, 2].axis('off')
        axarr[1, 2].axis('equal')

        fig.subplots_adjust(wspace=0.0)

        return fig, axarr

    def plot_test_C(self, gtype = 'vor', size = (10, 8), print_errors = True):
        """
        Analytical function comparison for the Helmholtz-Hodge decomposition of a vector field into
        its divergence-free and curl-free components.

        As the 'vor' mesh represents an open/free boundary, this condition best
        replicates the analytical math equations.

        """

        ax = self.centroid[0]
        ay = self.centroid[1]

        if gtype=='tri':
            xo = self.tri_verts[:, 0]
            yo = self.tri_verts[:, 1]

        elif gtype=='vor':
            xo = self.vor_verts[:,0]
            yo = self.vor_verts[:,1]

        else:
            raise Exception("'tri' and 'vor' are only valid mesh types.")

        # Define an analytical function representing a scalar potential:
        self.Phi_o = (xo - ax) ** 2 + (yo - ay) ** 2
        self.gPhix_o = 2 * (xo - ax) # analytical gradient of the scalar potential
        self.gPhiy_o = 2 * (yo - ay)

        # Define an analytical function representing a vector potential (oriented in the z-axis direction):
        self.Psi_o = xo * yo
        self.gPsix_o = xo # analytical gradient of the vector potential
        self.gPsiy_o = yo

        self.Fxo = self.gPhix_o - self.gPsiy_o # Construct an analytic vector field F = grad Phi + curl Psi
        self.Fyo = self.gPhiy_o + self.gPsix_o

        # Testing DEC-based Helmholtz-Hodge decomposition algorithm
        # Map the field to the mid points from the verts
        Fxom = self.verts_to_mids(self.Fxo, gtype=gtype)
        Fyom = self.verts_to_mids(self.Fyo, gtype=gtype)

        # Use the Helmholtz Hodge (HH) decomposition of the field on the DEC grids:
        cPsixm, cPsiym, gPhixm, gPhiym = self.helmholtz_hodge(Fxom, Fyom, gtype=gtype)

        self.gPhix_ = self.mids_to_verts(gPhixm, gtype=gtype)
        self.gPhiy_ = self.mids_to_verts(gPhiym, gtype=gtype)

        self.cPsix_ = self.mids_to_verts(cPsixm, gtype=gtype)
        self.cPsiy_ = self.mids_to_verts(cPsiym, gtype=gtype)

        # Attempt to reconstruct the force-field F using the DEC-derived HH components:
        self.Fx_ = self.gPhix_ + self.cPsix_
        self.Fy_ = self.gPhiy_ + self.cPsiy_


        self.error_Fx = np.sqrt((self.Fxo - self.Fx_)**2)
        self.error_Fy = np.sqrt((self.Fyo - self.Fy_)**2)

        self.test_errors = {'Fx': self.error_Fx.mean()/np.abs(self.Fxo).max(),
                            'Fy': self.error_Fy.mean()/np.abs(self.Fyo).max(),
                            }

        if print_errors:
            print("Percent error on estimates:")
            for k, v in self.test_errors.items():
                print(k, np.round(v*100, 4))


        # Plot the results:
        fig, axarr = plt.subplots(2, 3, figsize=size)

        axarr[0, 0].xaxis.set_ticklabels([])
        axarr[0, 0].yaxis.set_ticklabels([])

        axarr[0, 1].xaxis.set_ticklabels([])
        axarr[0, 1].yaxis.set_ticklabels([])

        axarr[0, 2].xaxis.set_ticklabels([])
        axarr[0, 2].yaxis.set_ticklabels([])

        axarr[1, 0].xaxis.set_ticklabels([])
        axarr[1, 0].yaxis.set_ticklabels([])

        axarr[1, 1].xaxis.set_ticklabels([])
        axarr[1, 1].yaxis.set_ticklabels([])

        axarr[1, 2].xaxis.set_ticklabels([])
        axarr[1, 2].yaxis.set_ticklabels([])

        axarr[0, 0].quiver(xo, yo, self.Fxo, self.Fyo, color='purple', zorder=10, headwidth = 5)
        axarr[0, 0].set_title('Original F (exact)')
        # axarr[0, 0].axis('tight')
        axarr[0, 0].axis('off')
        axarr[0, 0].axis('equal')
        # axarr[0, 0].axis(self.xyaxis)

        axarr[0, 1].quiver(xo, yo, self.gPhix_o, self.gPhiy_o, color='red', zorder=10, headwidth = 5)
        axarr[0, 1].set_title('Curl-free F (exact)')
        # axarr[0, 1].axis('tight')
        axarr[0, 1].axis('off')
        axarr[0, 1].axis('equal')
        # axarr[0, 1].axis(self.xyaxis)

        axarr[0, 2].quiver(xo, yo, -self.gPsiy_o, self.gPsix_o, color='blue', zorder=10, headwidth = 5)
        axarr[0, 2].set_title('Div-free F (exact)')
        # axarr[0, 2].axis('tight')
        axarr[0, 2].axis('off')
        axarr[0, 2].axis('equal')
        # axarr[0, 2].axis(self.xyaxis)

        axarr[1, 0].quiver(xo, yo, self.Fx_, self.Fy_, color='purple', zorder=10, headwidth = 5)
        axarr[1, 0].set_title('Reconstructed F (DEC)')
        # axarr[1, 0].axis('tight')
        axarr[1, 0].axis('off')
        axarr[1, 0].axis('equal')
        # axarr[1, 0].axis(self.xyaxis)

        axarr[1, 1].quiver(xo, yo, self.gPhix_, self.gPhiy_, color='red', zorder=10, headwidth = 5)
        axarr[1, 1].set_title('Curl-free F (DEC)')
        # axarr[1, 1].axis('tight')
        axarr[1, 1].axis('off')
        axarr[1, 1].axis('equal')
        # axarr[1, 1].axis(self.xyaxis)

        axarr[1, 2].quiver(xo, yo,  self.cPsix_, self.cPsiy_,  color='blue', zorder=10, headwidth = 5)
        axarr[1, 2].set_title('Div-free F (DEC)')
        # axarr[1, 2].axis('tight')
        axarr[1, 2].axis('off')
        axarr[1, 2].axis('equal')
        # axarr[1, 2].axis(self.xyaxis)

        fig.subplots_adjust(wspace=0.0)

    #
        return fig, axarr








