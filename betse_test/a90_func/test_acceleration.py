#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by equivalent influencers.
# See "LICENSE" for further details.

"""
Unit tests for the "KSL Sovereign" hardware acceleration module.
"""

import pytest
import numpy as np
from unittest.mock import patch
import sys
from betse.science.math import finitediff

# ....................{ TESTS                              }....................
def test_metal_laplacian_natural():
    """
    Scenario A: The "Happy Path".
    Runs the laplacian with whatever hardware is naturally available.
    """
    grid_size = 200
    F = np.random.rand(grid_size, grid_size).astype(np.float32)
    delx = 1.0

    # This should never crash, regardless of whether Torch is installed or not.
    result = finitediff.laplacian(F, delx)

    assert result.shape == F.shape
    assert not np.isnan(result).any()
    assert np.any(result != 0)


def test_forced_cpu_fallback():
    """
    Scenario B: The "Safety Simulation".
    Ensures that even if Torch IS installed, we can force the CPU path.
    If Torch is NOT installed, we verify the fallback works naturally.
    """
    grid_size = 200
    F = np.random.rand(grid_size, grid_size).astype(np.float32)
    delx = 1.0

    # Check if 'torch' is even installed in this test environment
    try:
        import torch
        # CASE 1: Torch exists. We must use Mocking to simulate a non-Metal machine.
        # We lie to the system and say "Metal is not available."
        with patch('torch.backends.mps.is_available', return_value=False):
            result_cpu = finitediff.laplacian(F, delx)
            
    except ImportError:
        # CASE 2: Torch is missing (Standard CI / MSI Laptop).
        # We don't need to mock anything. The code should naturally catch the
        # ImportError internally and use the Numpy fallback.
        result_cpu = finitediff.laplacian(F, delx)

    # Validation: The result must be valid
    assert result_cpu.shape == F.shape
    assert not np.isnan(result_cpu).any()
    assert np.any(result_cpu != 0)

