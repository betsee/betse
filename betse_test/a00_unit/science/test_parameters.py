#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :mod:`betse.science.parameters` submodule.
'''

# ....................{ IMPORTS                            }....................
from pytest import raises

# ....................{ TESTS                              }....................
def test_parameters_conf_file_missing(
    betse_temp_dir: 'py._path.local.LocalPath') -> None:
    '''
    Unit test that loading a simulation configuration whose top-level YAML file
    was copied *without* the external files it references raises a
    human-readable exception naming the missing file, the setting referencing
    that file, and the configuration file defining that setting *and* advising
    the end user on how to resolve this.

    This is the scenario reported by issue #6, in which an end user copied only
    the ``sample_sim/sample_sim.yaml`` file to another directory.

    Parameters
    ----------
    betse_temp_dir : py._path.local.LocalPath
        Object encapsulating a temporary directory isolated to this test.
    '''

    # Defer heavyweight imports.
    import shutil
    from betse.exceptions import BetseFileException
    from betse.science.parameters import Parameters

    # Absolute filename of a copy of the default simulation configuration file
    # in this temporary directory, copied *WITHOUT* the "extra_configs/"
    # subdirectory and other external files referenced by that file.
    conf_filepath = betse_temp_dir.join('Defiant.yaml')
    conf_filename = str(conf_filepath)
    shutil.copyfile(Parameters.conf_default_filename, conf_filename)

    # Assert that loading this configuration raises the expected exception.
    with raises(BetseFileException) as exception_info:
        Parameters.make(conf_filename=conf_filename)
    exception_message = str(exception_info.value)

    # Assert this exception message to name the missing file, the setting
    # referencing that file, and the configuration file defining that setting.
    assert 'expression_data.yaml' in exception_message
    assert '"general network/expression data file"' in exception_message
    assert conf_filename in exception_message

    # Assert this exception message to advise how to resolve this.
    assert 'cp -R' in exception_message
    assert 'betse config' in exception_message
