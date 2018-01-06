#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Create and electrodiffuses a suite of customizable general molecule in the BETSE ecosystem,
including functionality to pump the molecule, use it as a gating ligand, produce and consume it,
and use it an enzyme to facilitate another reaction. The molecule is assumed to be at low
concentrations and to not have a significant effect on system voltages or currents. This
module creates a structure containing all user-defined molecules, along with the facilities
to initialize, define the core computations for a simulation loop, save and report on data,
and plot.

"""

from betse.exceptions import BetsePyDotException
from betse.lib import libs


def set_net_opts(self, net_plot_opts, p):

    if net_plot_opts is not None:

        self.plot_net = net_plot_opts['plot network']
        self.net_font_name = net_plot_opts['font name']
        self.node_font_size = net_plot_opts['node font size']
        self.tit_font_size = net_plot_opts['sub-net font size']
        self.net_layout = net_plot_opts['layout']
        self.edge_width = net_plot_opts['edge linewidth']

        self.conc_node_font_color = net_plot_opts['substances']['font color']
        self.conc_node_color = net_plot_opts['substances']['node color']
        self.conc_shape = net_plot_opts['substances']['node shape']

        self.react_node_font_color = net_plot_opts['reactions']['font color']
        self.react_node_color = net_plot_opts['reactions']['node color']
        self.reaction_shape = net_plot_opts['reactions']['node shape']

        self.transp_node_font_color = net_plot_opts['transporters']['font color']
        self.transp_node_color = net_plot_opts['transporters']['node color']
        self.transporter_shape = net_plot_opts['transporters']['node shape']

        self.chan_node_font_color = net_plot_opts['channels']['font color']
        self.chan_node_color = net_plot_opts['channels']['node color']
        self.channel_shape = net_plot_opts['channels']['node shape']

        self.ed_node_font_color = net_plot_opts['electrodiffusion']['font color']
        self.ed_node_color = net_plot_opts['electrodiffusion']['node color']
        self.ed_shape = net_plot_opts['electrodiffusion']['node shape']

        self.vmem_node_font_color = net_plot_opts['Vmem']['font color']
        self.vmem_node_color = net_plot_opts['Vmem']['node color']
        self.vmem_shape = net_plot_opts['Vmem']['node shape']

        self.subnets_dicts = net_plot_opts['sub-nets']
        self.additional_edges = net_plot_opts.get('additional relationships', None)

    else: # set to defaults

        self.plot_net = p.plot_network
        self.net_font_name = 'DejaVu Sans'
        self.node_font_size = 16
        self.tit_font_size = 24
        self.net_layout = 'TB'
        self.edge_width = 2.0

        self.conc_node_font_color = 'Black'
        self.conc_node_color = 'LightCyan'
        self.conc_shape = 'ellipse'

        self.react_node_font_color = 'White'
        self.react_node_color = 'DarkSeaGreen'
        self.reaction_shape = 'rect'

        self.transp_node_font_color = 'White'
        self.transp_node_color = 'DarkSeaGreen'
        self.transporter_shape = 'diamond'

        self.chan_node_font_color = 'White'
        self.chan_node_color = 'DarkSeaGreen'
        self.channel_shape = 'pentagon'

        self.ed_node_font_color = 'White'
        self.ed_node_color = 'DarkSeaGreen'
        self.ed_shape = 'hexagon'

        self.vmem_node_font_color = 'White'
        self.vmem_node_color = 'maroon'
        self.vmem_shape = 'oval'

        self.subnets_dicts = None
        self.additional_edges = None

def plot_master_network(self, p):

    # If PyDot is unimportable, raise an exception.
    libs.die_unless_runtime_optional('pydot')

    # reserve import of pydot in case the user doesn't have it and needs to turn this functionality off:
    import pydot

    # create a graph object
    base_graph = pydot.Dot(graph_type='digraph', concentrate=False, nodesep=0.1, ranksep=0.3,
                                splines=True, strict = True, rankdir = self.net_layout)

    # add Vmem to the graph
    nde = pydot.Node('Vmem', style='filled', shape=self.vmem_shape, color = self.vmem_node_color,
                             fontcolor = self.vmem_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
    base_graph.add_node(nde)

    # add each substance as a node in the graph:
    for i, (name, val) in enumerate(self.cell_concs.items()):

        if name not in p.ions_dict:

            mol = self.molecules[name]

            nde = pydot.Node(name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
            base_graph.add_node(nde)

            if mol.simple_growth:
                # add node & edge for growth reaction component:
                rea_name = name + '_growth'
                rea_node = pydot.Node(rea_name, style='filled', color = self.react_node_color, shape=self.reaction_shape,
                             fontcolor = self.react_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                base_graph.add_node(rea_node)

                # if the substance has autocatalytic growth capacity add the edge in:
                base_graph.add_edge(pydot.Edge(rea_name, name, arrowhead='normal', coeff=1.0,
                                               penwidth = self.edge_width))

                # add node & edge for decay reaction component:
                rea_name = name + '_decay'

                rea_node = pydot.Node(rea_name, style='filled', color = self.react_node_color, shape=self.reaction_shape,
                             fontcolor = self.react_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                base_graph.add_node(rea_node)

                # if the substance has autocatalytic decay capacity add the edge in:
                base_graph.add_edge(pydot.Edge(name, rea_name, arrowhead='normal', coeff=1.0,
                                               penwidth = self.edge_width))

            if mol.ion_channel_gating:

                # if this substance gates for ion channels:

                # define a node corresponding to the ion channel:
                gated_node = pydot.Node(mol.gating_channel_name, style='filled',
                             color = self.chan_node_color, shape=self.channel_shape,
                             fontcolor = self.chan_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                base_graph.add_node(gated_node)

                # add the nodes for substance gating channel, if needed:
                if mol.gating_extracell is True:
                    substance_name = name + "_env"

                    nde = pydot.Node(name, style='filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                    base_graph.add_node(nde)

                else:
                    substance_name = name


                base_graph.add_edge(pydot.Edge(substance_name, gated_node, arrowhead='dot', color='blue',
                                               penwidth = self.edge_width))

                for ion_name in mol.gating_ion_name:

                    ion_Pmem = 'Pmem_' + ion_name

                    # add the edges for channel effect on ion concentration:
                    base_graph.add_edge(pydot.Edge(gated_node, ion_Pmem,  arrowhead='dot', color='blue',
                                               penwidth = self.edge_width))


        else:  # add the ion as a node, but don't worry about growth/decay

            nde = pydot.Node(name, style='filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
            base_graph.add_node(nde)

        # add expression nodes for electrodiffusion/diffusion:
        # check the Dmem value of the substance:
        dmem_check = self.Dmem[name]

        if dmem_check != 0.0:  # FIXME here is one place where Vmem can activate/inhibit the reaction itself...

            react_label = 'Pmem_'+ name # here, "Pmem_" stands for membrane permeability
            nde = pydot.Node(react_label, style='filled', color = self.ed_node_color, shape=self.ed_shape,
                             fontcolor = self.ed_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)

            # add the environmental complement for this substance:
            name_env = name + '_env'
            nde = pydot.Node(name_env, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)

            base_graph.add_edge(pydot.Edge(name, react_label, arrowhead='normal', coeff=1.0,
                                           penwidth = self.edge_width))

            base_graph.add_edge(pydot.Edge(react_label, name_env, arrowhead='normal', coeff=1.0,
                                           penwidth=self.edge_width))


    # deal with activators and inhibitors for substance growth------------------------------------------------
    for i, name in enumerate(self.molecules):

        mol = self.molecules[name]
        # add regulatory as nodes in the graph:

        if mol.simple_growth is True:

            gname = name + '_growth'

            graph_influencers(self, base_graph, gname, mol.growth_activators_list,
                mol.growth_inhibitors_list, p, reaction_zone = 'cell',
                zone_tags_a = mol.growth_activators_zone,
                zone_tags_i = mol.growth_inhibitors_zone)

        if mol.ion_channel_gating is True:

            graph_influencers(self, base_graph, mol.gating_channel_name, mol.ion_activators_list,
                mol.ion_inhibitors_list, p, reaction_zone = 'cell',
                zone_tags_a = mol.ion_activators_zone,
                zone_tags_i = mol.ion_inhibitors_zone)

    #------------------------------------------------------------------------------------------------------------------
    # Basic Bioelectric relations

    # detail how Vmem is affected via membrane permeability relationships to core ions:

    base_graph.add_edge(pydot.Edge('Pmem_Na', 'Vmem', arrowhead='dot', color='blue', penwidth = self.edge_width))

    base_graph.add_edge(pydot.Edge('Pmem_K', 'Vmem', arrowhead='tee', color='red', penwidth = self.edge_width))

    base_graph.add_edge(pydot.Edge('K_env', 'Vmem', arrowhead='dot', color='blue', penwidth = self.edge_width))

    #-------------------------------------------------------------------------------------------------------------------

    # if there are any reactions in the cytosol, add them to the graph
    if len(self.reactions) > 0:

        for i, name in enumerate(self.reactions):

            nde = pydot.Node(name, style='filled', color = self.react_node_color, shape=self.reaction_shape,
                             fontcolor = self.react_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)

    if len(self.reactions_env) > 0:

        for i, name in enumerate(self.reactions_env):
            nde = pydot.Node(name, style='filled', color=self.react_node_color, shape=self.reaction_shape,
                             fontcolor=self.react_node_font_color, fontname=self.net_font_name,
                             fontsize=self.node_font_size)

            base_graph.add_node(nde)


    # if there are any reactions in the mitochondria, add them to the graph
    if len(self.reactions_mit) > 0:

        for i, name in enumerate(self.reactions_mit):

            nde = pydot.Node(name, style='filled', color = self.react_node_color, shape=self.reaction_shape,
                             fontcolor = self.react_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)


    # if there are any env reactions, plot their edges on the graph--------------------------------------------------
    if len(self.reactions_env) > 0:

        for i, name in enumerate(self.reactions_env):

            rea = self.reactions_env[name]

            for i, react_name in enumerate(rea.reactants_list):
                rea_coeff = rea.reactants_coeff[i]
                react_name_e = react_name + '_env'

                nde = pydot.Node(react_name_e, style='filled', color=self.conc_node_color, shape=self.conc_shape,
                                 fontcolor=self.conc_node_font_color, fontname=self.net_font_name,
                                 fontsize=self.node_font_size)

                base_graph.add_node(nde)

                base_graph.add_edge(pydot.Edge(react_name_e, name, arrowhead='normal', coeff=rea_coeff,
                                               penwidth = self.edge_width))

            for j, prod_name in enumerate(rea.products_list):
                prod_coeff = rea.products_coeff[j]
                prod_name_e = prod_name + '_env'

                nde = pydot.Node(prod_name_e, style='filled', color=self.conc_node_color, shape=self.conc_shape,
                                 fontcolor=self.conc_node_font_color, fontname=self.net_font_name,
                                 fontsize=self.node_font_size)

                base_graph.add_node(nde)

                base_graph.add_edge(pydot.Edge(name, prod_name_e, arrowhead='normal', coeff=prod_coeff,
                                               penwidth = self.edge_width))

            graph_influencers(self, base_graph, name, rea.reaction_activators_list,
                                   rea.reaction_inhibitors_list, p, reaction_zone=rea.reaction_zone,
                                   zone_tags_a=rea.reaction_activators_zone,
                                   zone_tags_i=rea.reaction_inhibitors_zone)

    # if there are any reactions, plot their edges on the graph--------------------------------------------------
    if len(self.reactions) > 0:

        for i, name in enumerate(self.reactions):

            rea = self.reactions[name]

            for i, react_name in enumerate(rea.reactants_list):
                rea_coeff = rea.reactants_coeff[i]
                base_graph.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff=rea_coeff,
                                               penwidth=self.edge_width))

            for j, prod_name in enumerate(rea.products_list):
                prod_coeff = rea.products_coeff[j]
                base_graph.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff=prod_coeff,
                                               penwidth=self.edge_width))

            graph_influencers(self, base_graph, name, rea.reaction_activators_list,
                              rea.reaction_inhibitors_list, p, reaction_zone=rea.reaction_zone,
                              zone_tags_a=rea.reaction_activators_zone,
                              zone_tags_i=rea.reaction_inhibitors_zone)

    # if there are any mitochondria zone reactions, plot their edges on the graph (and react/prod nodes):
    if len(self.reactions_mit) > 0:

        for i, name in enumerate(self.reactions_mit):

            rea = self.reactions_mit[name]

            for i, react_name in enumerate(rea.reactants_list):
                react_name += '_mit'

                nde = pydot.Node(react_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                base_graph.add_node(nde)

                rea_coeff = rea.reactants_coeff[i]
                base_graph.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff=rea_coeff,
                                               penwidth = self.edge_width))

            for j, prod_name in enumerate(rea.products_list):
                prod_name += '_mit'

                nde = pydot.Node(prod_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
                base_graph.add_node(nde)

                prod_coeff = rea.products_coeff[j]
                base_graph.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff=prod_coeff,
                                               penwidth = self.edge_width))

            graph_influencers(self, base_graph, name, rea.reaction_activators_list,
                                   rea.reaction_inhibitors_list, p, reaction_zone=rea.reaction_zone,
                                   zone_tags_a=rea.reaction_activators_zone,
                                   zone_tags_i=rea.reaction_inhibitors_zone)

    # if there are any transporters, plot them on the graph:
    if len(self.transporters) > 0:

        for name in self.transporters:

            nde = pydot.Node(name, style='filled', color = self.transp_node_color, shape=self.transporter_shape,
                             fontcolor = self.transp_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)

        for name in self.transporters:

            trans = self.transporters[name]

            for i, (react_name, tag) in enumerate(zip(trans.reactants_list, trans.react_transport_tag)):

                rea_coeff = trans.reactants_coeff[i]

                if tag == 'cell_concs' or tag == 'mem_concs':

                    base_graph.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff=rea_coeff,
                                                   penwidth = self.edge_width))

                else:

                    if tag == 'env_concs':

                        react_name += '_env'

                    elif tag == 'mit_concs':

                        react_name += '_mit'

                    nde = pydot.Node(react_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                    base_graph.add_node(nde)

                    base_graph.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff=rea_coeff,
                                                   penwidth = self.edge_width))

            for j, (prod_name, tag) in enumerate(zip(trans.products_list, trans.prod_transport_tag)):

                prod_coeff = trans.products_coeff[j]

                if tag == 'cell_concs' or tag == 'mem_concs':

                    base_graph.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff=prod_coeff,
                                                   penwidth = self.edge_width))

                else:

                    if tag == 'env_concs':

                        prod_name += '_env'

                    elif tag == 'mit_concs':

                        prod_name += '_mit'

                    nde = pydot.Node(prod_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

                    base_graph.add_node(nde)

                    base_graph.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff=prod_coeff,
                                                   penwidth = self.edge_width))


            graph_influencers(self, base_graph, name, trans.transporter_activators_list,
                              trans.transporter_inhibitors_list, p, reaction_zone=trans.reaction_zone,
                              zone_tags_a=trans.transporter_activators_zone,
                              zone_tags_i=trans.transporter_inhibitors_zone)

    # if there are any channels, plot their type, ion  and Vmem relationships in the graph: -----------------------
    if len(self.channels) > 0:

        for i, name in enumerate(self.channels):

            chan = self.channels[name]

            chan_class = self.channels[name].channel_class
            channel_name = self.channels[name].channel_type

            if chan_class == 'Na':
                ion_name = ['Na']

            elif chan_class == 'K':
                ion_name = ['K']

            elif chan_class == 'Kir':
                ion_name = ['K']

            elif chan_class == 'Fun':
                ion_name = ['Na', 'K']

            elif chan_class == 'Ca':
                ion_name = ['Ca']

            elif chan_class == 'NaP':
                ion_name = ['Na']

            elif chan_class == 'Cat':
                ion_name = ['Na', 'K']

            elif chan_class == 'Cl':
                ion_name = ['Cl']

            else:
                ion_name = []

            # add the channel to the diagram:
            nde = pydot.Node(name, style='filled', color = self.chan_node_color, shape=self.channel_shape,
                             fontcolor = self.chan_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)

            base_graph.add_node(nde)

            for ion_n in ion_name:

                ion_Pmem = 'Pmem_' + ion_n

                base_graph.add_edge(pydot.Edge(name, ion_Pmem,  arrowhead='dot', color='blue',
                                               penwidth = self.edge_width))

            graph_influencers(self, base_graph, name, chan.channel_activators_list,
                                   chan.channel_inhibitors_list, p, reaction_zone='cell',
                                   zone_tags_a=chan.channel_activators_zone,
                                   zone_tags_i=chan.channel_inhibitors_zone)

    # if there are any user-defined relations, add them in:

    if self.additional_edges is not None:
        for node_list in self.additional_edges:

            node1, node2, relation = node_list

            if not base_graph.get_node(node1):

                raise BetsePyDotException('Additional node "{}" not found.'.format(node1))

            if not base_graph.get_node(node2):

                raise BetsePyDotException('Additional node "{}" not found.'.format(node2))

            if relation == 'activation':
                arrowh = 'dot'
                arrowc = 'blue'

            elif relation == 'inhibition':
                arrowh = 'tee'
                arrowc = 'red'

            else:
                arrowh = 'normal'
                arrowc = 'black'


            base_graph.add_edge(pydot.Edge(node1, node2, arrowhead=arrowh, color=arrowc,
                                           penwidth=self.edge_width))


    # if sub-graphs are defined, re-jigger the network to have clusters/sub-networks:
    if self.subnets_dicts is not None and self.subnets_dicts != 'None':

        if len(self.subnets_dicts) > 0:

            base_graph = make_subgraphs(self, base_graph, p)

    return base_graph

def make_subgraphs(self, base_graph, p):

    # If PyDot or Networkx are unimportable, raise an exception.
    libs.die_unless_runtime_optional('pydot', 'networkx')

    # reserve import of pydot and networkx in case the user doesn't have it and needs to turn this functionality off:
    import pydot
    import networkx as nx

    # convert the pydot-version base_graph into a networkx function so we can re-jigger nodes and make subgraphs...
    net_all = nx.from_pydot(base_graph)

    subnetworks = {}

    master_graph = pydot.Dot(graph_type='digraph', concentrate=False, nodesep=0.1, ranksep=0.3,
                                splines=True, strict = True, rankdir = self.net_layout)

    # Begin by writing all nodes and edges from the main network to the new pydot network
    master_graph = write_nodes(net_all, master_graph)
    master_graph = write_edges(net_all, master_graph)


    for subnet_opts in self.subnets_dicts:

        subnet_name = subnet_opts['name']
        subnet_nodes = subnet_opts['nodes']

        # expand subnetwork nodes list with nearest neighbours from main graph:
        extended_node_list = []

        for nde in subnet_nodes:

            if nde in net_all.nodes():
                extended_node_list.append(nde)
                nn = net_all.neighbors(nde)

                for n in nn:
                    extended_node_list.append(n)

        # reassign subgraph nodes to the extended list:
        subnet_nodes = extended_node_list

        tit_font = subnet_opts['title font color']
        box_color = subnet_opts['box shading color']

        subnet = nx.subgraph(net_all, subnet_nodes)
        subnetworks[subnet_name] = {'name': subnet_name, 'nx subnet': subnet, 'tit_font': tit_font,
                                    'box_color': box_color,  'nodes': subnet_nodes}



    for i, sn_dic in enumerate(subnetworks):

        cname = 'cluster_' + str(i)

        # declare subgraphs in pydot format:
        py_sub = pydot.Subgraph(cname, label=subnetworks[sn_dic]['name'], style='filled',
                                fontcolor = subnetworks[sn_dic]['tit_font'],
                            fontname=self.net_font_name, fontsize=self.tit_font_size,
                                color=subnetworks[sn_dic]['box_color'])

        subnetworks[sn_dic]['py subnet'] = py_sub

    # now each subnetwork has a networkx and empty pydot subgraph definition.
    # write each networkx subnet to the pydot subgraph:
    for i, sn_dic in enumerate(subnetworks):

        sn_py = write_nodes(subnetworks[sn_dic]['nx subnet'], subnetworks[sn_dic]['py subnet'])

        subnetworks[sn_dic]['py subnet'] = sn_py

        master_graph.add_subgraph(sn_py)

        master_graph = write_edges(subnetworks[sn_dic]['nx subnet'], master_graph)


    return master_graph

def graph_influencers(self, base_graph, name, a_list, i_list, p, reaction_zone = 'cell', zone_tags_a = None,
    zone_tags_i = None):
    """
    A helper method allowing activators and inhibitors, which have various options, to be added to a graph.

    Parameters
    --------------
    graph:    The pydot graph object to add to
    name             Name of the main molecule being activated or inhibited
    a_list:   List of activator names
    i_list:   List of inhibitor names
    reaction_zone:  Zone for the main reaction
    zone_tags_a:     list of zones activator works from
    zone_tags_i:     list of zones inhibitor work from.

    Returns
    --------
    graph           Updated pydot graph object

    """
    # If PyDot is unimportable, raise an exception.
    libs.die_unless_runtime_optional('pydot')

    import pydot

    if a_list != 'None' and a_list is not None:

        for act_name, zone_a in zip(a_list, zone_tags_a):

            independent_action = False

            if act_name.endswith('!'):  # if activator indicated to be independent
                # make the name in accordance with its actual identifier
                act_name = act_name[:-1]

            if reaction_zone == 'mit':

                act_name += '_mit'
                nde = pydot.Node(act_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
                base_graph.add_node(nde)

            if zone_a == 'env':

                act_name = act_name + '_env'
                nde = pydot.Node(act_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
                base_graph.add_node(nde)


            base_graph.add_edge(pydot.Edge(act_name, name, penwidth = self.edge_width, arrowhead='dot', color='blue'))


    if i_list != 'None' and i_list is not None:

        for inh_name, zone_i in zip(i_list, zone_tags_i):

            if reaction_zone == 'mit':

                inh_name += '_mit'
                nde = pydot.Node(inh_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
                base_graph.add_node(nde)

            if zone_i == 'env':

                inh_name = inh_name + '_env'
                nde = pydot.Node(inh_name, style = 'filled', color = self.conc_node_color, shape=self.conc_shape,
                             fontcolor = self.conc_node_font_color, fontname = self.net_font_name,
                             fontsize = self.node_font_size)
                base_graph.add_node(nde)

            base_graph.add_edge(
                pydot.Edge(inh_name, name, arrowhead='tee', color='red', penwidth = self.edge_width))

def write_nodes(nx_graph, py_graph):

    # If PyDot or Networkx are unimportable, raise an exception.
    libs.die_unless_runtime_optional('pydot')
    # reserve import of pydot and networkx in case the user doesn't have it and needs to turn this functionality off:
    import pydot

    for nx_node_name in nx_graph.node:

        nx_dic = nx_graph.node[nx_node_name]

        clr = nx_dic.get('color', None)
        fntclr = nx_dic.get('fontcolor', None)
        fntnm = nx_dic.get('fontname', None)
        fntsz = nx_dic.get('fontsize', None)
        shp = nx_dic.get('shape', None)
        stl = nx_dic.get('style', None)

        # if (clr is not None and fntclr is not None and fntnm is not None and shp is not None and
        #             shp is not None and stl is not None):

        nde = pydot.Node(nx_node_name, color=clr, fontcolor=fntclr, fontname=fntnm, fontsize=fntsz,
                         shape=shp, style=stl)

        py_graph.add_node(nde)

        # else:
        #
        #     logs.log_warning("WARNING! You have requested a sub-graph node that does not exist!\n"
        #                      " Ignoring request!")

    return py_graph

def write_edges(nx_graph, py_graph):

    # If PyDot or Networkx are unimportable, raise an exception.
    libs.die_unless_runtime_optional('pydot', 'networkx')
    # reserve import of pydot and networkx in case the user doesn't have it and needs to turn this functionality off:
    import pydot
    import networkx


    for na, nb in nx_graph.edges():

        edge_dic = nx_graph.edge[na][nb]

        arrow = edge_dic.get('arrowhead', None)
        linew = edge_dic.get('penwidth', None)
        clr = edge_dic.get('color', 'black')

        # if arrow is not None and linew is not None and clr is not None:

        eg = pydot.Edge(na, nb, arrowhead=arrow, penwidth=linew, color=clr)
        py_graph.add_edge(eg)

        # else:
        #
        #     logs.log_warning("WARNING! You have requested a sub-graph edge that does not exist!\n"
        #                      " Ignoring request!")

    return py_graph





