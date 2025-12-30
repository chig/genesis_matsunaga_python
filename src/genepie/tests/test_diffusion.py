# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
from .conftest import MSD_DATA
from .. import genesis_exe
import numpy as np


def test_diffusion_analysis():
    msd = np.loadtxt(MSD_DATA, dtype=np.float64)
    assert msd.ndim == 2, "MSD data should be 2D array"
    print(f"MSD data shape: {msd.shape}")

    ret = genesis_exe.diffusion_analysis(
            msd,
            time_step = 2,
            start = "20 %"
            )

    # Validate diffusion coefficient results
    assert ret is not None, "Diffusion result should not be None"
    assert len(ret) > 0, "Diffusion result should have at least one value"
    # Diffusion coefficients should be non-negative
    assert (ret >= 0).all(), "Diffusion coefficients should be non-negative"
    print(f"Diffusion coefficients (n={len(ret)}): {ret}")


def test_diffusion_analysis_zerocopy():
    """Test zerocopy interface."""
    msd = np.loadtxt(MSD_DATA, dtype=np.float64)
    assert msd.ndim == 2, "MSD data should be 2D array"
    print(f"MSD data shape: {msd.shape}")

    ndata = msd.shape[0]
    start_step = int(ndata * 0.2)  # 20% start

    result = genesis_exe.diffusion_analysis_zerocopy(
            msd,
            time_step=2.0,
            start_step=start_step
            )

    assert result is not None, "Result should not be None"
    assert result.out_data is not None, "out_data should not be None"
    assert result.diffusion_coefficients is not None, "diffusion_coefficients should not be None"
    assert (result.diffusion_coefficients >= 0).all(), "Diffusion coefficients should be non-negative"
    print(f"Diffusion coefficients (zerocopy): {result.diffusion_coefficients}")


def test_diffusion_analysis_zerocopy_full():
    """Test zerocopy_full interface (pre-allocated result arrays)."""
    msd = np.loadtxt(MSD_DATA, dtype=np.float64)
    assert msd.ndim == 2, "MSD data should be 2D array"
    print(f"MSD data shape: {msd.shape}")

    ndata = msd.shape[0]
    start_step = int(ndata * 0.2)  # 20% start

    # Run zerocopy version for comparison
    result_zerocopy = genesis_exe.diffusion_analysis_zerocopy(
            msd,
            time_step=2.0,
            start_step=start_step
            )

    # Run zerocopy_full version
    result_full = genesis_exe.diffusion_analysis_zerocopy_full(
            msd,
            time_step=2.0,
            start_step=start_step
            )

    # Verify shapes match
    assert result_full.out_data.shape == result_zerocopy.out_data.shape, \
           f"out_data shape mismatch: {result_full.out_data.shape} vs {result_zerocopy.out_data.shape}"
    assert result_full.diffusion_coefficients.shape == result_zerocopy.diffusion_coefficients.shape, \
           f"diffusion_coefficients shape mismatch"

    # Verify values match exactly
    max_diff_out = np.max(np.abs(result_full.out_data - result_zerocopy.out_data))
    max_diff_coeff = np.max(np.abs(result_full.diffusion_coefficients - result_zerocopy.diffusion_coefficients))
    print(f"Max diff out_data: {max_diff_out}")
    print(f"Max diff diffusion_coefficients: {max_diff_coeff}")

    np.testing.assert_allclose(result_full.out_data, result_zerocopy.out_data,
                               rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result_full.diffusion_coefficients, result_zerocopy.diffusion_coefficients,
                               rtol=1e-10, atol=1e-10)
    print("zerocopy_full results match zerocopy results!")


def main():
    test_diffusion_analysis()
    test_diffusion_analysis_zerocopy()
    test_diffusion_analysis_zerocopy_full()


if __name__ == "__main__":
    main()
