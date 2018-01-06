#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.lib.yaml` subpackage.
'''

# ....................{ IMPORTS                            }....................

# ....................{ TESTS                              }....................
def test_yaml_copy(
    betse_sim_conf: 'SimConfTestInternal',
    betse_temp_dir: 'LocalPath',
) -> None:
    '''
    Test the capacity of the :mod:`betse.lib.yaml.yamls` submodule to reliably
    copy a deserialized simulation configuration (including both the top-level
    YAML file for this configuration *and* all external paths internally
    referenced and hence required by this file).

    Parameters
    ----------
    betse_sim_conf : SimConfTestInternal
        Object encapsulating a temporary simulation configuration file.
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current test.
    '''

    # Absolute path of a subdirectory with arbitrary basename residing in this
    # temporary directory.
    new_subdirpath = betse_temp_dir.join('Ten_Forward')

    # Create this subdirectory
    new_subdirpath.ensure(dir=True)

    # Absolute path of a new simulation configuration file with arbitrary
    # basename to be created in this subdirectory.
    new_sim_conf_filepath = new_subdirpath.join('El-Aurian_Guinan.yaml')

    # Absolute path of this file as a string rather than "LocalPath" object.
    new_sim_conf_filename = str(new_sim_conf_filepath)

    # Simulation configuration loaded from this file.
    p = betse_sim_conf.p

    # Copy this configuration to this subdirectory.
    p.save(new_sim_conf_filename)

    # Assert this file to have been created.
    assert new_sim_conf_filepath.check(file=1)


def test_yaml_roundtrip(betse_sim_conf: 'SimConfTestInternal') -> None:
    '''
    Test the capacity of the :mod:`betse.lib.yaml.yamls` submodule to reliably
    roundtrip (i.e., load and save without loss) the default YAML-formatted
    simulation configuration file, whose complexity warrants explicit testing.

    Parameters
    ----------
    betse_sim_conf : SimConfTestInternal
        Object encapsulating a temporary simulation configuration file.
    '''

    # Simulation configuration loaded from this file.
    p = betse_sim_conf.p

    # Absolute path of this file.
    p_conf_filename = p.conf_filename

    #FIXME: Test addition, deletion, and modification of non-scalar list items.

    # Modify an arbitrary boolean setting of this configuration.
    p.is_ecm = not p.is_ecm
    p_sim_ECM_expected = p.is_ecm

    # Modify an arbitrary numeric setting of this configuration.
    p.cell_polarizability *= 3.1415
    p_cell_polarizability_expected = p.cell_polarizability

    # Modify an arbitrary string setting of this configuration.
    p.seed_pickle_basename += '~'
    p_seed_pickle_basename_expected = p.seed_pickle_basename

    # Save these changes back to the same file.
    p.save_inplace()

    # Close this file.
    p.unload()

    # Read these changes back into the same in-memory object.
    p.load(p_conf_filename)

    # Ensure these changes were roundtripped across this I/O.
    assert p.is_ecm == p_sim_ECM_expected
    assert p.cell_polarizability == p_cell_polarizability_expected
    assert p.seed_pickle_basename == p_seed_pickle_basename_expected
