#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.lib.yaml` subpackage.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_yaml_roundtrip(betse_sim_config: 'SimTestState') -> None:
    '''
    Test the capacity of the :mod:`betse.lib.yaml.yamls` submodule to reliably
    roundtrip (i.e., load and save without loss) the default YAML-formatted
    simulation configuration file, whose complexity warrants explicit testing.

    Parameters
    ----------
    betse_sim_config : SimTestState
        Object encapsulating a temporary simulation configuration file.
    '''

    # Simulation configuration loaded from this file.
    p = betse_sim_config.p

    # Absolute path of this file.
    p_conf_filename = p.conf_filename

    #FIXME: Test addition, deletion, and modification of non-scalar list items.

    # Modify an arbitrary boolean setting of this configuration.
    p.sim_ECM = not p.sim_ECM
    p_sim_ECM_expected = p.sim_ECM

    # Modify an arbitrary numeric setting of this configuration.
    p.cell_polarizability *= 3.1415
    p_cell_polarizability_expected = p.cell_polarizability

    # Modify an arbitrary string setting of this configuration.
    p.pickle_seed_basename += '~'
    p_pickle_seed_basename_expected = p.pickle_seed_basename

    # Save these changes back to the same file.
    p.overwrite()

    # Close this file.
    p.unread()

    # Read these changes back into the same in-memory object.
    p.read(p_conf_filename)

    # Ensure these changes were roundtripped across this I/O.
    assert p.sim_ECM == p_sim_ECM_expected
    assert p.cell_polarizability == p_cell_polarizability_expected
    assert p.pickle_seed_basename == p_pickle_seed_basename_expected
